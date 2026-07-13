"""Stage 2 — Functional characterization of variants at three certainty tiers.

Data sources:
    MEASURED
      - MaveDB DMS scores (api.mavedb.org)
      - ClinVar functional records with PS3/BS3 codes

    PREDICTED
      - AlphaMissense per-variant class + score
      - ESM (skip in default run; heavy)
      - REVEL (skip unless dbNSFP provided)

    ANNOTATED
      - UniProt features (active/binding/disulfide/mutagen/...)
      - InterPro domains
      - AlphaFold pLDDT (per-residue confidence)
      - PDB interfaces (skip in default; needs PDB parser)

Every variant × source table is written to parquet with a `certainty_tier` column.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

log = logging.getLogger("rcvtc.stage2")


def run_stage2(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage2"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage2 start | %d genes", len(gene_ids))

    from .stage2a_mavedb import fetch_mavedb_scores
    from .stage2b_alphamissense import fetch_alphamissense_scores
    from .stage2c_uniprot import fetch_uniprot_features
    from .stage2d_alphafold import fetch_alphafold_plddt

    counts = {}

    if cfg["stage2"]["measured"].get("mavedb", {}).get("enabled"):
        log.info("2a: MaveDB DMS scores")
        df = fetch_mavedb_scores(cfg, gene_ids, outdir / "2a_mavedb")
        counts["mavedb_scores"] = len(df) if df is not None else 0

    if cfg["stage2"]["predicted"].get("alphamissense", {}).get("enabled"):
        log.info("2b: AlphaMissense scores")
        df = fetch_alphamissense_scores(cfg, gene_ids, outdir / "2b_alphamissense")
        counts["alphamissense_scores"] = len(df) if df is not None else 0

    if cfg["stage2"]["annotated"].get("uniprot_features", {}).get("enabled"):
        log.info("2c: UniProt features")
        df = fetch_uniprot_features(cfg, gene_ids, outdir / "2c_uniprot")
        counts["uniprot_features"] = len(df) if df is not None else 0

    if cfg["stage2"]["annotated"].get("alphafold_plddt", {}).get("enabled"):
        log.info("2d: AlphaFold pLDDT")
        df = fetch_alphafold_plddt(cfg, gene_ids, outdir / "2d_alphafold")
        counts["alphafold_residues"] = len(df) if df is not None else 0

    log.info("stage2 counts: %s", counts)
    return {"counts": counts}
