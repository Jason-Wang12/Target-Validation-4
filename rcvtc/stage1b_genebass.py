"""Stage 1b — Genebass burden scan across 4,131 UKB phenotype partitions.

Critical implementation constraint: PyArrow dataset-level filter pushdown across
the full 4,131-partition Hive layout times out on the S3 FUSE mount. Instead, we
iterate partitions sequentially and filter each in-memory. Each partition is
~780 KB / ~15K rows and reads in ~1 second, so the full scan for a small gene
set completes in ~1-2 hours per annotation mask.

The output is:
- one row per (gene, annotation, phenotype) combination that meets the p-value threshold
- with tier tag = INFERRED (statistical test on GB summary stats)
"""
from __future__ import annotations
import logging
import os
import time
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage1b")

BASE_LOF = "/mnt/datalake/genebass/genebass_pLoF_all_filtered.parquet"
BASE_MIS = "/mnt/datalake/genebass/genebass_missense_LC_all_filtered.parquet"


def scan_genebass_burden(cfg: dict, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    """Sequential per-partition scan filtered to config genes."""
    outdir.mkdir(parents=True, exist_ok=True)
    stage1 = cfg["stage1"]
    if not stage1.get("associations", {}).get("genebass", {}).get("enabled", False):
        log.info("genebass disabled in config")
        return None

    annotations = stage1["associations"]["genebass"].get("annotations", ["pLoF", "missense|LC"])
    gene_set = set(gene_ids.keys())
    log.info("genebass scan | genes=%s | annotations=%s", sorted(gene_set), annotations)

    # Suggestive threshold for early filter; final calls apply Bonferroni downstream
    suggestive_p = 1.0e-4

    all_rows: list[pd.DataFrame] = []

    for base, tag in [(BASE_LOF, "pLoF"), (BASE_MIS, "missense|LC")]:
        if tag not in annotations:
            continue
        if not os.path.isdir(base):
            log.error("genebass parquet root missing: %s", base)
            continue
        partitions = sorted(
            d for d in os.listdir(base)
            if d.startswith("phenotype_slug=")
        )
        log.info("genebass %s: %d partitions to scan", tag, len(partitions))

        t0 = time.time()
        kept = 0
        skipped_read_err = 0
        for i, part in enumerate(partitions):
            slug = part.replace("phenotype_slug=", "")
            path = f"{base}/{part}/part-00000.parquet"
            try:
                df = pd.read_parquet(path)
            except Exception as e:
                skipped_read_err += 1
                if skipped_read_err <= 5:
                    log.warning("read failed for %s: %s", part, e)
                continue
            sub = df[df["gene"].isin(gene_set) & (df["Pvalue"] <= suggestive_p)]
            if len(sub) > 0:
                sub = sub.copy()
                sub["phenotype_slug"] = slug
                sub["annotation_mask"] = tag
                all_rows.append(sub)
                kept += len(sub)

            # Progress log fires unconditionally (was previously nested under 'kept' -> silent)
            if (i + 1) % 500 == 0:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                remain = (len(partitions) - i - 1) / rate
                log.info(
                    "  %s: %d/%d partitions (kept=%d) | %.1f/s | ~%.0fs left",
                    tag, i + 1, len(partitions), kept, rate, remain,
                )
        log.info(
            "  %s complete | kept=%d | read_errors=%d | wall=%.1fs",
            tag, kept, skipped_read_err, time.time() - t0,
        )

    if not all_rows:
        log.warning("genebass: no suggestive hits for any gene")
        empty = pd.DataFrame(columns=[
            "gene", "annotation_mask", "phenotype_slug", "pheno_description",
            "BETA_Burden", "SE_Burden", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
            "certainty_tier", "source",
        ])
        safe_write_parquet(empty, outdir / "genebass_associations.parquet")
        return empty

    out = pd.concat(all_rows, ignore_index=True)
    out = out.rename(columns={"annotation": "annotation_col"})   # avoid schema collision
    # tag every row for downstream tier-audit
    out["certainty_tier"] = Tier.INFERRED.value
    out["source"] = SOURCES.get("genebass", "genebass_ukb_wes_450k")
    out = out.sort_values(["gene", "annotation_mask", "Pvalue"]).reset_index(drop=True)
    safe_write_parquet(out, outdir / "genebass_associations.parquet")
    log.info("genebass_associations.parquet written | rows=%d", len(out))
    return out
