"""Stage 1a — gnomAD v4 variant catalog per gene.

Uses the gnomAD GraphQL API to enumerate coding variants in each gene's
canonical transcript. Applies:
    - popmax AF < config.stage1.af_filter.max_popmax_af
    - filter flags NOT in config.stage1.af_filter.exclude_flags
    - min AC / min coverage as configured

All gnomAD allele frequencies are MEASURED (empirical count in cohort).
"""
from __future__ import annotations
import logging
import time
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage1a")

URL = "https://gnomad.broadinstitute.org/api"

# Note: we ask for `variants` at gene level with dataset gnomad_r4 (v4).
# The GraphQL schema has been stable; if this breaks we surface it clearly.
GQL = """
query GeneVariants($sym: String!, $ref: ReferenceGenomeId!, $ds: DatasetId!) {
  gene(gene_symbol: $sym, reference_genome: $ref) {
    gene_id
    symbol
    canonical_transcript_id
    variants(dataset: $ds) {
      variant_id
      rsid
      pos
      ref
      alt
      consequence
      hgvsc
      hgvsp
      flags
      exome {
        ac
        an
        af
        filters
        populations { id ac an }
      }
      genome {
        ac
        an
        af
        filters
        populations { id ac an }
      }
    }
  }
}
"""

# Consequences we consider coding for pipeline scope (LoF + missense).
CODING = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost",
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
}


def _popmax_af(pops: list[dict] | None) -> float:
    """Return max AF across non-ancestry-super populations."""
    if not pops:
        return 0.0
    # gnomAD 'id' values that are super-populations (afr, amr, asj, eas, fin,
    # mid, nfe, sas). Exclude the pseudo-populations and sex-strata.
    keep = {"afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"}
    afs = []
    for p in pops:
        pid = (p.get("id") or "").lower()
        if pid not in keep:
            continue
        ac, an = p.get("ac"), p.get("an")
        if ac is None or an is None or an == 0:
            continue
        afs.append(ac / an)
    return max(afs) if afs else 0.0


def fetch_gnomad_variants(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    if not cfg["stage1"].get("variant_source", "").startswith("gnomad"):
        log.info("stage1.variant_source is not gnomad; skipping")
        return None

    af_cfg = cfg["stage1"]["af_filter"]
    max_af = float(af_cfg.get("max_popmax_af", 1e-3))
    exclude_flags = set(af_cfg.get("exclude_flags", ["AC0", "InbreedingCoeff", "LowQual", "RF"]))

    sess = get_session()
    out_rows = []
    for sym, g in gene_ids.items():
        t0 = time.time()
        r = sess.post(URL, json={
            "query": GQL,
            "variables": {"sym": sym, "ref": "GRCh38", "ds": "gnomad_r4"},
        }, timeout=90)
        if not r.ok:
            log.warning("gnomAD failed | %s | %s", sym, r.status_code)
            continue
        js = r.json()
        if "errors" in js:
            log.warning("gnomAD errors | %s | %s", sym, js["errors"][:2])
        data = (js.get("data") or {}).get("gene")
        if not data:
            log.warning("gnomAD: no gene data for %s", sym)
            continue
        variants = data.get("variants") or []
        log.info("gnomAD | %s | %d raw variants | %.1fs", sym, len(variants), time.time() - t0)

        kept = 0
        for v in variants:
            conseq = v.get("consequence")
            if conseq not in CODING:
                continue

            ex = v.get("exome") or {}
            ge = v.get("genome") or {}

            # Filter flags from either callset
            flags_ex = set(ex.get("filters") or [])
            flags_ge = set(ge.get("filters") or [])
            if (flags_ex & exclude_flags) or (flags_ge & exclude_flags):
                continue
            variant_flags = set(v.get("flags") or [])

            popmax = max(_popmax_af(ex.get("populations")), _popmax_af(ge.get("populations")))
            # Total AC / AN
            ac_total = (ex.get("ac") or 0) + (ge.get("ac") or 0)
            an_total = (ex.get("an") or 0) + (ge.get("an") or 0)
            af_total = (ac_total / an_total) if an_total else 0.0

            if ac_total < af_cfg.get("min_popmax_ac", 1):
                continue
            if popmax >= max_af:
                continue

            out_rows.append({
                "gene_symbol": sym,
                "ensembl_id": g.ensembl_id,
                "variant_id": v["variant_id"],
                "rsid": v.get("rsid"),
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "consequence": conseq,
                "hgvsc": v.get("hgvsc"),
                "hgvsp": v.get("hgvsp"),
                "ac_total": ac_total,
                "an_total": an_total,
                "af_total": af_total,
                "af_popmax": popmax,
                "variant_flags": ";".join(sorted(variant_flags)) if variant_flags else None,
                "exome_filters": ";".join(sorted(flags_ex)) if flags_ex else None,
                "genome_filters": ";".join(sorted(flags_ge)) if flags_ge else None,
                "certainty_tier": Tier.MEASURED.value,
                "source": SOURCES.get("gnomad", "gnomad_v4_1"),
            })
            kept += 1

        log.info("  %s: %d kept (rare coding, passing filters)", sym, kept)

    if not out_rows:
        log.warning("gnomAD: no variants collected")
        empty = pd.DataFrame(columns=[
            "gene_symbol", "ensembl_id", "variant_id", "rsid", "pos", "ref", "alt",
            "consequence", "hgvsc", "hgvsp", "ac_total", "an_total", "af_total", "af_popmax",
            "variant_flags", "exome_filters", "genome_filters", "certainty_tier", "source",
        ])
        safe_write_parquet(empty, outdir / "gnomad_variants.parquet")
        return empty

    df = pd.DataFrame(out_rows).sort_values(["gene_symbol", "consequence", "af_popmax"], ascending=[True, True, False]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "gnomad_variants.parquet")
    log.info("gnomad_variants.parquet written | rows=%d", len(df))
    return df
