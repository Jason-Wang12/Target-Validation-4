"""Stage 5 — FinnGen R12 replication (INFERRED tier).

For each gene in the run, we pull the top single-variant associations at the
gene locus from FinnGen's public GraphQL-free JSON API:

    /api/gene_phenos/{gene}            — top variant per phenotype at locus
    /api/gene_functional_variants/{gene} — functional missense/LoF variants
                                            with all significant phenotypes

These are **single-variant** effects — not the same test as Genebass gene-burden
— so we treat the concordance as replication when the direction & magnitude
suggest the same biology, not exact p-value matching.

Downstream (stage 9) uses this to flag whether a gene-phenotype pair from
Genebass has any FinnGen support.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage5")

BASE = "https://r12.finngen.fi/api"


def fetch_finngen(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()
    if not cfg.get("stage5", {}).get("finngen_r12", {}).get("enabled", True):
        log.info("FinnGen R12 disabled")
        return pd.DataFrame()

    all_phenos = []
    all_variants = []

    for sym in gene_ids.keys():
        # ── gene_phenos: top-per-phenotype hits at the gene locus ──────────
        r = sess.get(f"{BASE}/gene_phenos/{sym}", timeout=45)
        if not r.ok:
            log.warning("FinnGen gene_phenos failed | %s | %s", sym, r.status_code)
        else:
            payload = r.json() or {}
            for ph in payload.get("phenotypes", []):
                a = ph.get("assoc", {}) or {}
                v = ph.get("variant", {}) or {}
                all_phenos.append({
                    "gene_symbol": sym,
                    "phenocode": a.get("phenocode"),
                    "phenostring": a.get("phenostring"),
                    "category": a.get("category"),
                    "beta": a.get("beta"),
                    "sebeta": a.get("sebeta"),
                    "pval": a.get("pval"),
                    "mlogp": a.get("mlogp"),
                    "maf": a.get("maf"),
                    "maf_case": a.get("maf_case"),
                    "maf_control": a.get("maf_control"),
                    "n_case": a.get("n_case"),
                    "n_control": a.get("n_control"),
                    "n_sample": a.get("n_sample"),
                    "variant_chrom": v.get("chr") or v.get("chrom"),
                    "variant_pos": v.get("pos"),
                    "variant_ref": v.get("ref"),
                    "variant_alt": v.get("alt"),
                    "variant_varid": v.get("varid"),
                    "variant_annotation": v.get("annotation"),
                    "certainty_tier": Tier.INFERRED.value,
                    "source": SOURCES.get("finngen", "finngen_r12"),
                })

        # ── gene_functional_variants: missense/LoF variants + phenos ──────
        r2 = sess.get(f"{BASE}/gene_functional_variants/{sym}", timeout=45)
        if not r2.ok:
            log.warning("FinnGen gene_functional_variants failed | %s | %s", sym, r2.status_code)
        else:
            payload2 = r2.json() or []
            for row in payload2:
                var = row.get("var", {}) or {}
                for ph in row.get("significant_phenos", []) or []:
                    all_variants.append({
                        "gene_symbol": sym,
                        "rsids": row.get("rsids"),
                        "variant_chrom": var.get("chr") or var.get("chrom"),
                        "variant_pos": var.get("pos"),
                        "variant_ref": var.get("ref"),
                        "variant_alt": var.get("alt"),
                        "annotation": var.get("most_severe") or var.get("annotation"),
                        "gene_most_severe": var.get("gene_most_severe"),
                        "info": var.get("info"),
                        "phenocode": ph.get("phenocode"),
                        "phenostring": ph.get("phenostring"),
                        "category": ph.get("category"),
                        "beta": ph.get("beta"),
                        "sebeta": ph.get("sebeta"),
                        "pval": ph.get("pval"),
                        "mlogp": ph.get("mlogp"),
                        "maf": ph.get("maf"),
                        "maf_case": ph.get("maf_case"),
                        "maf_control": ph.get("maf_control"),
                        "n_case": ph.get("n_case"),
                        "n_control": ph.get("n_control"),
                        "n_sample": ph.get("n_sample"),
                        "certainty_tier": Tier.INFERRED.value,
                        "source": SOURCES.get("finngen", "finngen_r12"),
                    })

        log.info(
            "FinnGen | %s | %d gene_phenos rows, %d functional-variant rows",
            sym,
            sum(1 for p in all_phenos if p["gene_symbol"] == sym),
            sum(1 for p in all_variants if p["gene_symbol"] == sym),
        )

    phenos_df = pd.DataFrame(all_phenos)
    var_df = pd.DataFrame(all_variants)

    # Coerce variant_annotation and info fields to string (mixed dict/str/None)
    for df_ in (phenos_df, var_df):
        for col in ("variant_annotation", "annotation", "info", "gene_most_severe"):
            if col in df_.columns:
                df_[col] = df_[col].astype(str).where(df_[col].notna(), None)

    safe_write_parquet(phenos_df, outdir / "finngen_gene_phenos.parquet")
    safe_write_parquet(var_df, outdir / "finngen_functional_variants.parquet")
    log.info("finngen: %d gene_phenos rows, %d variant rows",
             len(phenos_df), len(var_df))
    return phenos_df


def run_stage5(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage5"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage5 start | %d genes", len(gene_ids))

    df = fetch_finngen(cfg, gene_ids, outdir / "5a_finngen")
    counts = {"finngen_gene_phenos": len(df) if df is not None else 0}

    log.info("stage5 counts: %s", counts)
    return {"counts": counts}
