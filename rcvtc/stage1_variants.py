"""Stage 1 — variant catalog (gnomAD v4) + rare-variant PheWAS/burden associations.

Sub-stages:
    1a. gnomAD v4 variant catalog per gene (canonical transcript, rare filter).
    1b. Genebass burden scan across 4,131 phenotypes (per-partition read).
    1c. ExPheWas gene-level PheWAS via v1 REST.
    1d. Open Targets `evidences(datasourceIds:["gene_burden"])` pass-through.

Every emitted association is tagged with a certainty tier. Genebass/ExPheWas/OT
associations are INFERRED (statistical test on summary stats). gnomAD population
allele frequencies are MEASURED (empirical count in reference cohort).
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Any

from .certainty import Tier, Claim, ClaimSet, SOURCES
from .utils.logging import stage_logger

log = stage_logger("stage1")


def run_stage1(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage1"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    genes = list(gene_ids.keys())
    log.info("stage1 start | %d genes | %s", len(genes), genes)

    # Delegate to sub-stage modules (imported lazily so we can build each in isolation)
    from .stage1a_gnomad import fetch_gnomad_variants
    from .stage1b_genebass import scan_genebass_burden
    from .stage1c_expheway import fetch_expheway_phewas
    from .stage1d_ot_burden import fetch_ot_gene_burden

    counts: dict[str, int] = {}

    log.info("1a: gnomAD variant catalog")
    variants_df = fetch_gnomad_variants(cfg, gene_ids, outdir / "1a_variants")
    counts["gnomad_variants"] = len(variants_df) if variants_df is not None else 0

    log.info("1b: Genebass burden scan")
    genebass_df = scan_genebass_burden(cfg, gene_ids, outdir / "1b_genebass")
    counts["genebass_associations"] = len(genebass_df) if genebass_df is not None else 0

    log.info("1c: ExPheWas gene-level PheWAS")
    expheway_df = fetch_expheway_phewas(cfg, gene_ids, outdir / "1c_expheway")
    counts["expheway_associations"] = len(expheway_df) if expheway_df is not None else 0

    log.info("1d: Open Targets gene_burden")
    ot_df = fetch_ot_gene_burden(cfg, gene_ids, outdir / "1d_ot_burden")
    counts["ot_gene_burden_evidences"] = len(ot_df) if ot_df is not None else 0

    log.info("stage1 counts: %s", counts)
    return {"counts": counts}
