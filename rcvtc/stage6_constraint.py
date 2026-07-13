"""Stage 6 — Constraint & mechanism verdict.

Sources of constraint (all MEASURED tier — count-based observations from gnomAD):
    - LOEUF = oe_lof_upper (upper 90 % CI of observed/expected LoF)
    - obs_lof, exp_lof, obs_mis, exp_mis, obs_syn, exp_syn
    - pLI, lof_z, mis_z, syn_z
    - oe_lof, oe_mis, oe_syn (with 90 % CI where available)

Mechanism verdict rule (per PLAN.md v3, stage6.gof_signal_rules):

  GoF suggested when ALL of:
    - Any Genebass phenotype with |missense|LC BETA| >= |pLoF BETA|
      (from stage 3's mask_comparison output)
    - LOEUF > 0.6 (gene not haploinsufficient overall)
    - oe_mis < 0.6 (missense strongly selected against — activating variants tolerated)
    - Existing curated GoF term in stage 4 curated evidence (optional)

  LoF haploinsufficient when:
    - Genebass pLoF BETA magnitude > missense|LC BETA magnitude (concordant)
    - LOEUF < 0.6 (haploinsufficient) OR any pathogenic LoF ClinVar variant
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage6")

GNOMAD_GQL = "https://gnomad.broadinstitute.org/api"

CONSTRAINT_QUERY = """
query C($sym: String!) {
  gene(gene_symbol: $sym, reference_genome: GRCh38) {
    symbol
    gnomad_constraint {
      exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper
      exp_mis obs_mis oe_mis oe_mis_lower oe_mis_upper
      exp_syn obs_syn oe_syn
      pLI lof_z mis_z syn_z
    }
  }
}
"""


def fetch_constraint(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()

    rows = []
    for sym in gene_ids.keys():
        r = sess.post(GNOMAD_GQL,
                      json={"query": CONSTRAINT_QUERY, "variables": {"sym": sym}},
                      timeout=30)
        if not r.ok:
            log.warning("gnomAD constraint failed | %s | %s", sym, r.status_code)
            continue
        d = r.json()
        c = ((d.get("data") or {}).get("gene") or {}).get("gnomad_constraint") or {}
        if not c:
            log.warning("no constraint data for %s", sym)
            continue

        rows.append({
            "gene_symbol": sym,
            "loeuf": c.get("oe_lof_upper"),
            "oe_lof": c.get("oe_lof"),
            "oe_lof_lower": c.get("oe_lof_lower"),
            "obs_lof": c.get("obs_lof"),
            "exp_lof": c.get("exp_lof"),
            "oe_mis": c.get("oe_mis"),
            "oe_mis_lower": c.get("oe_mis_lower"),
            "oe_mis_upper": c.get("oe_mis_upper"),
            "obs_mis": c.get("obs_mis"),
            "exp_mis": c.get("exp_mis"),
            "oe_syn": c.get("oe_syn"),
            "pLI": c.get("pLI"),
            "lof_z": c.get("lof_z"),
            "mis_z": c.get("mis_z"),
            "syn_z": c.get("syn_z"),
            "certainty_tier": Tier.MEASURED.value,  # count-based observations from population WES
            "source": SOURCES.get("gnomad", "gnomad_v4_1"),
        })
        log.info(
            "constraint | %s | LOEUF=%.2f | pLI=%.2g | oe_mis=%.2f",
            sym, c.get("oe_lof_upper") or float('nan'),
            c.get("pLI") or float('nan'),
            c.get("oe_mis") or float('nan'),
        )

    if not rows:
        empty = pd.DataFrame()
        safe_write_parquet(empty, outdir / "gene_constraint.parquet")
        return empty

    df = pd.DataFrame(rows)
    safe_write_parquet(df, outdir / "gene_constraint.parquet")
    return df


import re

_LOF_TOKENS = ("*", "ter", "fs", "del", "dup")


def _classify_protein_change(pc: str) -> str:
    """Classify a single protein_change string into missense / lof / unknown."""
    if pc is None or pd.isna(pc):
        return "unknown"
    s = str(pc).lower()
    tok = s.split(",")[0].strip()
    if any(t in tok for t in _LOF_TOKENS):
        return "lof"
    if re.match(r"[a-z]+\d+[a-z*]", tok):
        return "missense"
    return "unknown"


def _mechanism_verdict(sym: str, con_row: pd.Series | None,
                       mask_rows: pd.DataFrame | None,
                       clinvar: pd.DataFrame | None,
                       curated: pd.DataFrame | None,
                       cfg: dict) -> dict:
    """Produce a single-line verdict dict for one gene.

    Priority of signals:
      1. Stage 3 mask patterns from Genebass (best evidence when present)
      2. ClinVar pathogenic protein-change breakdown (missense-only vs mixed vs LoF-heavy)
      3. gnomAD constraint (LOEUF/oe_mis)
      4. Curated evidence GoF terms
    """
    rules = cfg.get("stage6", {}).get("gof_signal_rules", {}) or {}
    loeuf_min_gof = rules.get("constraint_pattern_gof", {}).get("loeuf_min", 0.6)
    oe_mis_max_gof = rules.get("constraint_pattern_gof", {}).get("missense_oe_max", 0.6)

    # Extract constraint numbers
    loeuf = con_row["loeuf"] if con_row is not None else None
    oe_mis = con_row["oe_mis"] if con_row is not None else None
    mis_z = con_row["mis_z"] if con_row is not None and "mis_z" in con_row else None

    # Extract mask signal
    gof_mask_signal = False
    lof_mask_signal = False
    strongest_mask_pattern = None
    if mask_rows is not None and not mask_rows.empty:
        g = mask_rows[mask_rows["gene"] == sym] if "gene" in mask_rows.columns else mask_rows.copy()
        if not g.empty and "mask_pattern" in g.columns:
            counts = g["mask_pattern"].value_counts()
            gof_mask_signal = counts.get("gof_activating", 0) > 0
            lof_mask_signal = counts.get("lof_haploinsufficiency", 0) > 0
            strongest_mask_pattern = counts.idxmax() if len(counts) else None

    # ClinVar pathogenic protein-change breakdown
    n_path_missense = 0
    n_path_lof = 0
    if clinvar is not None and not clinvar.empty:
        cv_g = clinvar[(clinvar["gene_symbol"] == sym) & (clinvar["clinsig_bucket"] == "pathogenic_lp")]
        if not cv_g.empty:
            classes = cv_g["protein_change"].map(_classify_protein_change)
            n_path_missense = int((classes == "missense").sum())
            n_path_lof = int((classes == "lof").sum())
    clinvar_path_lof = n_path_lof > 0
    # Missense-only pattern: any pathogenic missense, zero pathogenic LoF, and >=1 evidence
    # Even 1 missense with strong tier-2 review + GoF corroboration is meaningful (e.g. JAK2 V617F).
    clinvar_missense_only = (n_path_missense >= 1) and (n_path_lof == 0)

    # Curated GoF term detection
    gof_curation_term = False
    if curated is not None and not curated.empty and "gene_symbol" in curated.columns:
        cur_g = curated[curated["gene_symbol"] == sym]
        # combine any text-y curated fields for pattern matching
        text_cols = [c for c in ["clinical_significances", "variant_functional_consequence_label",
                                 "disease_name", "variant_aa_descriptions"] if c in cur_g.columns]
        curated_text = " ".join(
            " ".join(str(v) for v in cur_g[c].dropna().astype(str).tolist()) for c in text_cols
        ).lower()
        gof_curation_term = any(t in curated_text for t in
                                ["gain of function", "activating mutation", "constitutively active",
                                 "gain_of_function"])

    # === Verdict logic (evidence-weighted) ===
    #
    # We combine multiple orthogonal signals rather than short-circuiting on
    # the first-available one. Burden mask patterns are informative but can
    # mislead for GoF genes (e.g. JAK2): a single driver missense (V617F) gets
    # diluted by non-driver missense carriers in the burden average, so the
    # missense β appears smaller than the pLoF β even though the missense
    # signal is biologically dominant. Therefore mask-derived LoF signal is
    # NOT decisive against a strong ClinVar/curated GoF signature.
    verdict = "undetermined"
    reasons = []

    # STRONG GoF override: missense-only ClinVar pathogenic + GoF curation term
    # or LOEUF-tolerant gene. This beats any mask signal (see docstring above).
    strong_gof_override = (
        clinvar_missense_only and
        (gof_curation_term or (loeuf is not None and loeuf >= loeuf_min_gof))
    )
    if strong_gof_override:
        verdict = "gof_activating"
        reasons.append(f"pathogenic ClinVar variants all missense ({n_path_missense}), no LoF")
        if gof_curation_term:
            reasons.append("GoF term present in curated evidence")
        if loeuf is not None and loeuf >= loeuf_min_gof:
            reasons.append(f"LOEUF={loeuf:.2f} (gene tolerant to LoF, consistent with GoF/DN mechanism)")
        if lof_mask_signal:
            reasons.append("Genebass burden pattern showed pLoF β > missense|LC β — likely burden-averaging artifact (single missense driver diluted by non-driver carriers)")

    # STRONG LoF: unambiguous mask signal + ClinVar LoF alleles present
    elif lof_mask_signal and clinvar_path_lof and n_path_lof >= 3:
        verdict = "lof_haploinsufficiency"
        reasons.append("Genebass pLoF BETA > missense|LC BETA (mask cross-comparison)")
        reasons.append(f"ClinVar pathogenic variants include {n_path_lof} LoF alleles")
        if loeuf is not None:
            reasons.append(f"LOEUF={loeuf:.2f}")

    # Mask signal only (no ClinVar cross-check)
    elif gof_mask_signal and not lof_mask_signal:
        verdict = "gof_activating"
        reasons.append("Genebass missense|LC BETA >= pLoF BETA (mask cross-comparison)")
    elif lof_mask_signal and not gof_mask_signal:
        verdict = "lof_haploinsufficiency"
        reasons.append("Genebass pLoF BETA > missense|LC BETA (mask cross-comparison)")

    # ClinVar-only fallbacks
    if verdict == "undetermined":
        if clinvar_missense_only:
            verdict = "missense_dominant"
            reasons.append(f"all pathogenic missense ({n_path_missense}) — GoF or dominant-negative candidate")
        elif clinvar_path_lof and n_path_lof >= max(3, n_path_missense * 0.2):
            verdict = "lof_haploinsufficiency"
            reasons.append(f"ClinVar pathogenic variants include {n_path_lof} LoF alleles")
            if loeuf is not None and loeuf < loeuf_min_gof:
                reasons.append(f"LOEUF={loeuf:.2f} confirms haploinsufficiency")

    # Constraint-only fallback
    if verdict == "undetermined" and loeuf is not None:
        if loeuf < 0.35:
            verdict = "lof_haploinsufficiency"
            reasons.append(f"LOEUF={loeuf:.2f} (strong LoF constraint)")
        elif loeuf >= loeuf_min_gof and oe_mis is not None and oe_mis <= oe_mis_max_gof:
            verdict = "gof_candidate"
            reasons.append(f"constraint pattern GoF-compatible (LOEUF={loeuf:.2f}, oe_mis={oe_mis:.2f})")

    if verdict == "undetermined" and clinvar_path_lof:
        verdict = "lof_probable"
        reasons.append("ClinVar pathogenic LoF variants present (single-source)")

    return {
        "gene_symbol": sym,
        "loeuf": loeuf,
        "oe_mis": oe_mis,
        "mis_z": mis_z,
        "gof_mask_signal": gof_mask_signal,
        "lof_mask_signal": lof_mask_signal,
        "strongest_mask_pattern": strongest_mask_pattern,
        "n_path_missense": n_path_missense,
        "n_path_lof": n_path_lof,
        "gof_curation_term": gof_curation_term,
        "clinvar_path_lof": clinvar_path_lof,
        "clinvar_missense_only": clinvar_missense_only,
        "verdict": verdict,
        "reasons": " | ".join(reasons) if reasons else None,
        "certainty_tier": Tier.INFERRED.value,
        "source": "rcvtc_mechanism_verdict_v1",
    }


def run_stage6(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage6"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage6 start | %d genes", len(gene_ids))

    # 6a: pull constraint metrics
    con_df = fetch_constraint(cfg, gene_ids, outdir / "6a_constraint")

    # Load inputs from earlier stages
    mask_path = Path(cfg["outputs_dir"]) / "stage3" / "mask_comparison.parquet"
    clinvar_path = Path(cfg["outputs_dir"]) / "stage4" / "4a_clinvar" / "clinvar_variants.parquet"
    curated_path = Path(cfg["outputs_dir"]) / "stage4" / "4b_ot_curated" / "ot_curated_evidences.parquet"

    mask_rows = pd.read_parquet(mask_path) if mask_path.exists() else None
    clinvar = pd.read_parquet(clinvar_path) if clinvar_path.exists() else None
    curated = pd.read_parquet(curated_path) if curated_path.exists() else None

    verdicts = []
    for sym in gene_ids.keys():
        cr = con_df[con_df["gene_symbol"] == sym].iloc[0] if not con_df.empty and (con_df["gene_symbol"] == sym).any() else None
        v = _mechanism_verdict(sym, cr, mask_rows, clinvar, curated, cfg)
        verdicts.append(v)
        log.info("verdict | %s | %s", sym, v["verdict"])

    vdf = pd.DataFrame(verdicts)
    safe_write_parquet(vdf, outdir / "mechanism_verdicts.parquet")
    log.info("mechanism_verdicts.parquet written | rows=%d", len(vdf))

    return {"counts": {"gene_constraint": len(con_df), "verdicts": len(vdf)}}
