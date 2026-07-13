"""Stage 1c — ExPheWas gene-level PheWAS via v1 REST API.

Root: https://exphewas.statgen.org/v1/api
Pattern (verified in preflight):
    GET /gene/name/{symbol}     -> {ensembl_id, name, has_results, ...}
    GET /gene/{ENSG}/results    -> [{outcome_id, outcome_label, p, nlog10p, bonf, ...}]

All ExPheWas rows are INFERRED (gene-level logistic/linear PheWAS on UKB).
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

log = logging.getLogger("rcvtc.stage1c")

BASE = "https://exphewas.statgen.org/v1/api"


def fetch_expheway_phewas(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    stage1 = cfg["stage1"]
    if not stage1.get("associations", {}).get("expheway", {}).get("enabled", False):
        log.info("expheway disabled in config")
        return None

    sess = get_session()
    subset = stage1["associations"]["expheway"].get("analysis_subset", "BOTH")
    analysis = stage1["associations"]["expheway"].get("analysis_type", "PHECODES")
    timeout_s = stage1["associations"]["expheway"].get("timeout_s_per_query", 30)

    all_rows = []
    for sym, g in gene_ids.items():
        ensg = g.ensembl_id
        # Gene metadata + has_results check (must use /gene/ensembl/<id>, not /gene/<id>)
        r_meta = sess.get(f"{BASE}/gene/ensembl/{ensg}", timeout=timeout_s)
        if not r_meta.ok:
            log.warning("expheway meta failed | %s | %s", sym, r_meta.status_code)
            continue
        meta = r_meta.json()
        if not meta.get("has_results", False):
            log.warning("expheway has no results for %s", sym)
            continue

        # Results
        t0 = time.time()
        r = sess.get(f"{BASE}/gene/{ensg}/results", timeout=timeout_s)
        if not r.ok:
            log.warning("expheway results failed | %s | %s", sym, r.status_code)
            continue
        rows = r.json()
        log.info("expheway | %s | %d rows | %.1fs", sym, len(rows), time.time() - t0)
        for row in rows:
            if row.get("analysis_subset") != subset:
                continue
            if row.get("analysis_type") != analysis:
                continue
            row["gene_symbol"] = sym
            row["ensembl_id"] = ensg
            row["certainty_tier"] = Tier.INFERRED.value
            row["source"] = SOURCES.get("expheway", "expheway_v1")
            all_rows.append(row)

    if not all_rows:
        log.warning("expheway: no rows collected")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows).sort_values(["gene_symbol", "p"]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "expheway_associations.parquet")
    log.info("expheway_associations.parquet written | rows=%d", len(df))
    return df
