"""Stage 2c — UniProt features (ANNOTATED tier).

Feature types requested (configurable):
    ACT_SITE, BINDING, SITE, LIPID, MOD_RES, DISULFID, MUTAGEN, VARIANT

Note: UniProt's `type` values in JSON output are `Active site`, `Binding site`,
`Site`, `Lipidation`, `Modified residue`, `Disulfide bond`, `Mutagenesis`,
`Natural variant` — we handle this normalization in the module.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage2c")

UNIPROT = "https://rest.uniprot.org/uniprotkb"

# Map from the short codes (UniProt flat-file style) to JSON `type` strings
CODE_TO_JSON_TYPE = {
    "ACT_SITE":  "Active site",
    "BINDING":   "Binding site",
    "SITE":      "Site",
    "LIPID":     "Lipidation",
    "MOD_RES":   "Modified residue",
    "DISULFID":  "Disulfide bond",
    "MUTAGEN":   "Mutagenesis",
    "VARIANT":   "Natural variant",
    "TRANSMEM":  "Transmembrane",
    "DOMAIN":    "Domain",
    "REGION":    "Region",
    "MOTIF":     "Motif",
}


def _flatten_location(loc: dict) -> tuple[int | None, int | None]:
    """Extract start, end positions from a UniProt location object."""
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    return start, end


def fetch_uniprot_features(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()
    codes = cfg["stage2"]["annotated"]["uniprot_features"].get(
        "feature_types",
        ["ACT_SITE", "BINDING", "SITE", "LIPID", "MOD_RES", "DISULFID", "MUTAGEN", "VARIANT"],
    )
    wanted_types = {CODE_TO_JSON_TYPE.get(c, c) for c in codes}
    log.info("UniProt features | wanted types: %s", sorted(wanted_types))

    all_rows = []
    for sym, g in gene_ids.items():
        if not g.uniprot_id:
            log.warning("skip %s: no UniProt ID", sym)
            continue
        r = sess.get(f"{UNIPROT}/{g.uniprot_id}.json", timeout=30)
        if not r.ok:
            log.warning("UniProt fetch failed | %s | %s", sym, r.status_code)
            continue
        d = r.json()
        feats = d.get("features", [])
        kept = 0
        for f in feats:
            if f.get("type") not in wanted_types:
                continue
            start, end = _flatten_location(f.get("location") or {})
            all_rows.append({
                "gene_symbol": sym,
                "uniprot_id": g.uniprot_id,
                "feature_type": f.get("type"),
                "description": f.get("description"),
                "start": start,
                "end": end,
                "evidences": ";".join(
                    e.get("evidenceCode", "") for e in (f.get("evidences") or [])
                ) or None,
                "certainty_tier": Tier.ANNOTATED.value,
                "source": SOURCES.get("uniprot", "uniprotkb"),
            })
            kept += 1
        log.info("  %s (%s): %d features kept", sym, g.uniprot_id, kept)

    if not all_rows:
        log.warning("UniProt: no features collected")
        safe_write_parquet(pd.DataFrame(), outdir / "uniprot_features.parquet")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows).sort_values(["gene_symbol", "start"]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "uniprot_features.parquet")
    log.info("uniprot_features.parquet written | rows=%d", len(df))
    return df
