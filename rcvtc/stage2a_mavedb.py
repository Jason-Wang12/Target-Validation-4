"""Stage 2a — MaveDB DMS scores (MEASURED tier).

Uses the correct API host `api.mavedb.org` (not `www.mavedb.org`). The search
uses ScoreSetsSearch schema: `{"targets": [gene_symbols], "published": True}`.

Score CSV columns (verified against urn:mavedb:00001269-a-1 = LDLR LDL uptake):
    accession, hgvs_nt, hgvs_splice, hgvs_pro, score, hgvsp, aapos, sd, se, df

The `score` field is study-normalized (assay-specific). Callers should NOT
compare scores across scoresets without renormalization.
"""
from __future__ import annotations
import io
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage2a")

BASE = "https://api.mavedb.org"


def fetch_mavedb_scores(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()
    all_rows = []
    scoreset_meta = []

    for sym, g in gene_ids.items():
        r = sess.post(f"{BASE}/api/v1/score-sets/search",
                      json={"targets": [sym], "published": True}, timeout=30)
        if not r.ok:
            log.warning("MaveDB search failed for %s | %s", sym, r.status_code)
            continue
        d = r.json()
        sets = d.get("scoreSets") if isinstance(d, dict) else (d if isinstance(d, list) else [])
        # Filter to those where the exact target name matches (case-insensitive)
        matched = [s for s in sets if any(
            (t.get("name", "").upper() == sym.upper())
            for t in s.get("targetGenes", [])
        )]
        log.info("MaveDB | %s | %d matched scoresets (of %d returned)",
                 sym, len(matched), len(sets))
        for s in matched:
            urn = s["urn"]
            # Pull scores
            r2 = sess.get(f"{BASE}/api/v1/score-sets/{urn}/scores", timeout=90)
            if not r2.ok:
                log.warning("MaveDB scores fetch failed | %s | %s", urn, r2.status_code)
                continue
            try:
                df = pd.read_csv(io.StringIO(r2.text))
            except Exception as e:
                log.warning("MaveDB parse failed | %s | %s", urn, e)
                continue
            df["gene_symbol"] = sym
            df["scoreset_urn"] = urn
            df["scoreset_title"] = s.get("title")
            df["certainty_tier"] = Tier.MEASURED.value
            df["source"] = SOURCES.get("mavedb", "mavedb")
            all_rows.append(df)
            scoreset_meta.append({
                "gene_symbol": sym,
                "scoreset_urn": urn,
                "title": s.get("title"),
                "published_date": s.get("publishedDate"),
                "num_variants": s.get("numVariants"),
                "short_description": s.get("shortDescription"),
            })
            log.info("  %s | %s | %d variant scores", sym, urn, len(df))

    if scoreset_meta:
        meta_df = pd.DataFrame(scoreset_meta)
        safe_write_parquet(meta_df, outdir / "mavedb_scoresets.parquet")

    if not all_rows:
        log.warning("MaveDB: no scores collected for any gene")
        safe_write_parquet(pd.DataFrame(), outdir / "mavedb_scores.parquet")
        return pd.DataFrame()

    out = pd.concat(all_rows, ignore_index=True)
    safe_write_parquet(out, outdir / "mavedb_scores.parquet")
    log.info("mavedb_scores.parquet written | rows=%d", len(out))
    return out
