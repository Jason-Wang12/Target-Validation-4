"""Stage 4 — Published allelic-series retrieval + curated evidence.

Data sources:
    ANNOTATED
      - ClinVar 2+ star pathogenic/benign variants (E-utils)
      - Open Targets curated evidence (eva, orphanet, genomics_england,
        gene2phenotype, uniprot_literature)

Literature searches are performed at synthesis time (Stage 9) by the driving
agent because they require the LiteratureSearch tool, not deterministic API calls.
"""
from __future__ import annotations
import logging
from pathlib import Path

log = logging.getLogger("rcvtc.stage4")


def run_stage4(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage4"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage4 start | %d genes", len(gene_ids))

    from .stage4a_clinvar import fetch_clinvar_variants
    from .stage4b_ot_curated import fetch_ot_curated

    counts = {}

    if cfg.get("stage4", {}).get("clinvar", {}).get("enabled"):
        log.info("4a: ClinVar variants (pathogenic + benign, star-filtered)")
        df = fetch_clinvar_variants(cfg, gene_ids, outdir / "4a_clinvar")
        counts["clinvar_variants"] = len(df) if df is not None else 0

    if cfg.get("stage4", {}).get("ot_curated", {}).get("enabled"):
        log.info("4b: Open Targets curated evidence")
        df = fetch_ot_curated(cfg, gene_ids, outdir / "4b_ot_curated")
        counts["ot_curated_evidences"] = len(df) if df is not None else 0

    log.info("stage4 counts: %s", counts)
    return {"counts": counts}
