"""Stage 2b — AlphaMissense per-variant scores via Ensembl VEP.

Rather than downloading the ~5GB AlphaMissense table, we call Ensembl's VEP REST
endpoint with `AlphaMissense=1`. This gives us am_class + am_pathogenicity per
variant per transcript.

Input: the stage 1a gnomAD variants parquet.
Output: variant_id -> (am_class, am_pathogenicity), canonical transcript only.

AlphaMissense scores are PREDICTED tier (computational model).
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

log = logging.getLogger("rcvtc.stage2b")

VEP_BASE = "https://rest.ensembl.org"


def fetch_alphamissense_scores(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    # Load stage 1a variants
    v_path = Path(cfg["outputs_dir"]) / "stage1" / "1a_variants" / "gnomad_variants.parquet"
    if not v_path.exists():
        log.warning("Stage 1a variants not found at %s; skipping AlphaMissense", v_path)
        return pd.DataFrame()
    variants = pd.read_parquet(v_path)
    if len(variants) == 0:
        log.warning("Stage 1a variants file is empty")
        return pd.DataFrame()

    # We only need missense_variant for AlphaMissense
    miss = variants[variants["consequence"] == "missense_variant"].copy()
    log.info("AlphaMissense | %d missense variants across %d genes",
             len(miss), miss.gene_symbol.nunique())
    if len(miss) == 0:
        return pd.DataFrame()

    # POST batches to Ensembl VEP with region-strings
    # Regions: "1 55039974 55039974 G/T" is the accepted format
    sess = get_session()
    all_rows = []
    batch_size = 200                     # VEP's documented per-request cap
    hgvs_list = miss["variant_id"].tolist()   # already in "chr-pos-ref-alt" form

    # Convert variant_id -> VEP HGVS-like input string (VEP accepts "chr pos . ref alt . . .")
    # Simpler: post as region strings
    regions = []
    for vid in hgvs_list:
        parts = vid.split("-")
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts
        regions.append(f"{chrom} {pos} {int(pos) + len(ref) - 1} {ref}/{alt} 1")
    log.info("posting %d regions to Ensembl VEP in batches of %d", len(regions), batch_size)

    # Map region-string -> gene_symbol/variant_id (order is preserved)
    reg_meta = list(zip(hgvs_list, miss["gene_symbol"].tolist()))

    for batch_start in range(0, len(regions), batch_size):
        batch = regions[batch_start:batch_start + batch_size]
        batch_meta = reg_meta[batch_start:batch_start + batch_size]

        t0 = time.time()
        r = sess.post(
            f"{VEP_BASE}/vep/human/region",
            json={"variants": batch, "AlphaMissense": 1, "canonical": 1},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=60,
        )
        if not r.ok:
            log.warning("VEP batch %d failed | %s | body=%s", batch_start, r.status_code, r.text[:400])
            continue
        results = r.json()
        # Match result back to variant by returned input string
        by_input = {res.get("input"): res for res in results}
        for reg, (vid, sym) in zip(batch, batch_meta):
            res = by_input.get(reg)
            if not res:
                continue
            for tc in res.get("transcript_consequences", []):
                if tc.get("gene_symbol") != sym:
                    continue
                if tc.get("canonical") != 1:
                    continue
                am = tc.get("alphamissense") or {}
                all_rows.append({
                    "gene_symbol": sym,
                    "variant_id": vid,
                    "transcript_id": tc.get("transcript_id"),
                    "amino_acids": tc.get("amino_acids"),
                    "protein_start": tc.get("protein_start"),
                    "am_class": am.get("am_class"),
                    "am_pathogenicity": am.get("am_pathogenicity"),
                    "certainty_tier": Tier.PREDICTED.value,
                    "source": SOURCES.get("alphamissense", "alphamissense_v1"),
                })
        log.info("  batch %d/%d done in %.1fs",
                 batch_start // batch_size + 1,
                 (len(regions) + batch_size - 1) // batch_size,
                 time.time() - t0)

    if not all_rows:
        log.warning("AlphaMissense: no variants annotated")
        safe_write_parquet(pd.DataFrame(), outdir / "alphamissense_scores.parquet")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows).sort_values(["gene_symbol", "protein_start"]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "alphamissense_scores.parquet")
    log.info("alphamissense_scores.parquet written | rows=%d", len(df))
    return df
