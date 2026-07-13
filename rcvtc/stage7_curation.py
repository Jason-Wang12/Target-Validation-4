"""Stage 7 — Curated evidence synthesis (genetics-first weighting).

Consumes stage 4 outputs (ClinVar + OT curated) and produces a per-gene evidence
ledger with:

  - Per-datasource evidence counts (ClinVar, EVA, Orphanet, Genomics England,
    Gene2Phenotype, UniProt literature)
  - Weighted per-disease evidence score (genetics-first weighting profile)
  - Top-ranked disease per gene
  - Optional literature-count enrichment (populated at synthesis time)

All rows are ANNOTATED tier — human-curated relationships aggregated with a
transparent numeric weighting scheme.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage7")

# Genetics-first weight profile (per plan v3)
# Higher weight = more decisive evidence.
DEFAULT_WEIGHTS = {
    "eva": 1.0,                 # ClinVar via EVA - direct variant assertion
    "orphanet": 0.9,            # gene-disease + inheritance
    "genomics_england": 0.85,   # PanelApp green/amber gene panels
    "gene2phenotype": 0.85,     # DDG2P inheritance-annotated
    "uniprot_literature": 0.7,  # UniProt-linked publications
    # ClinVar star2+ ann evidence gets its own scoring below
}

STAR_WEIGHT = {0: 0.2, 1: 0.5, 2: 1.0, 3: 1.4, 4: 1.6}


def _score_ot(row) -> float:
    """OT evidence weighted by datasource weight × OT reported score."""
    w = DEFAULT_WEIGHTS.get(row["datasource_id"], 0.5)
    s = row.get("score")
    return float(s) * w if pd.notna(s) else 0.0


def _score_clinvar(row) -> float:
    """ClinVar variant weighted by review-stars and pathogenic/benign bucket."""
    if row["clinsig_bucket"] != "pathogenic_lp":
        # benign counts less; still recorded but with reduced weight
        return -0.5 * STAR_WEIGHT.get(int(row["stars"]), 0.2)
    return 1.0 * STAR_WEIGHT.get(int(row["stars"]), 0.2)


def run_stage7(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage7"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage7 start | %d genes", len(gene_ids))

    ot_path = Path(cfg["outputs_dir"]) / "stage4" / "4b_ot_curated" / "ot_curated_evidences.parquet"
    cv_path = Path(cfg["outputs_dir"]) / "stage4" / "4a_clinvar" / "clinvar_variants.parquet"

    ot_df = pd.read_parquet(ot_path) if ot_path.exists() else pd.DataFrame()
    cv_df = pd.read_parquet(cv_path) if cv_path.exists() else pd.DataFrame()

    # Weighted evidence scores at the OT (gene, disease) level
    per_disease_rows = []
    if not ot_df.empty:
        ot_df = ot_df.copy()
        ot_df["weighted_score"] = ot_df.apply(_score_ot, axis=1)
        agg = ot_df.groupby(["gene_symbol", "disease_id", "disease_name"], dropna=False).agg(
            n_evidence=("evidence_id", "count"),
            max_score=("score", "max"),
            median_score=("score", "median"),
            weighted_sum=("weighted_score", "sum"),
            weighted_max=("weighted_score", "max"),
            datasources=("datasource_id", lambda x: ";".join(sorted(set(x)))),
        ).reset_index()
        agg["certainty_tier"] = Tier.ANNOTATED.value
        agg["source"] = "opentargets_curated_v25.06"
        per_disease_rows.append(agg)

    per_disease = pd.concat(per_disease_rows, ignore_index=True) if per_disease_rows else pd.DataFrame()
    safe_write_parquet(per_disease, outdir / "curated_per_disease.parquet")

    # Gene-level ledger — counts of pathogenic ClinVar star2+, curated datasources
    gene_rows = []
    for sym in gene_ids.keys():
        row = {"gene_symbol": sym}
        # ClinVar counts by stars/bucket
        if not cv_df.empty and (cv_df["gene_symbol"] == sym).any():
            g = cv_df[cv_df["gene_symbol"] == sym]
            path = g[g["clinsig_bucket"] == "pathogenic_lp"]
            benign = g[g["clinsig_bucket"] == "benign_lb"]
            row["clinvar_pathogenic_n"] = int(len(path))
            row["clinvar_benign_n"] = int(len(benign))
            row["clinvar_max_stars"] = int(g["stars"].max())
            row["clinvar_score"] = float(g.apply(_score_clinvar, axis=1).sum())
        # OT per-datasource counts
        if not ot_df.empty and (ot_df["gene_symbol"] == sym).any():
            og = ot_df[ot_df["gene_symbol"] == sym]
            for ds in DEFAULT_WEIGHTS.keys():
                row[f"ot_{ds}_n"] = int(len(og[og["datasource_id"] == ds]))
            row["ot_weighted_sum"] = float(og["weighted_score"].sum())
            row["ot_n_distinct_diseases"] = int(og["disease_id"].nunique())
        # Top-ranked disease per gene (by weighted_sum then weighted_max)
        if not per_disease.empty and (per_disease["gene_symbol"] == sym).any():
            g = per_disease[per_disease["gene_symbol"] == sym].sort_values(
                ["weighted_sum", "weighted_max"], ascending=False
            ).head(3)
            row["top_diseases"] = ";".join(
                f"{r['disease_name']}({r['weighted_sum']:.2f})" for _, r in g.iterrows()
            )
        row["certainty_tier"] = Tier.ANNOTATED.value
        row["source"] = "rcvtc_curated_synthesis_v1"
        gene_rows.append(row)

    gene_df = pd.DataFrame(gene_rows)
    safe_write_parquet(gene_df, outdir / "curated_gene_ledger.parquet")
    log.info(
        "stage7 wrote | per_disease=%d rows | gene_ledger=%d rows",
        len(per_disease), len(gene_df),
    )
    for _, r in gene_df.iterrows():
        log.info("ledger | %s | clinvar_path=%s | ot_weighted=%s | top=%s",
                 r["gene_symbol"],
                 r.get("clinvar_pathogenic_n", "-"),
                 (f"{r.get('ot_weighted_sum'):.1f}" if r.get("ot_weighted_sum") is not None else "-"),
                 r.get("top_diseases", "-"))

    return {"counts": {"curated_per_disease": len(per_disease), "curated_gene_ledger": len(gene_df)}}
