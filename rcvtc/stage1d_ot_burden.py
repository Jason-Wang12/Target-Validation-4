"""Stage 1d — Open Targets gene_burden evidence (AZ PheWAS 470K pass-through).

Open Targets ingests the AZ PheWAS gene-level burden results and exposes them
via the `gene_burden` datasource on evidence rows (pre-filtered to p < 1e-7 by
OT's own ingestion pipeline). Since AZ PheWAS does not offer a live JSON API,
this is the closest we can get to the AZ 470K burden numbers programmatically.

The OT schema requires BOTH `ensemblIds` and `efoIds` on the `evidences` field,
so the flow is:
    1. target(ensemblId) -> associatedDiseases (filtered to genetic_association) -> disease.id list
    2. For each disease, pull evidences(datasourceIds:["gene_burden"], ensemblIds, efoIds)

All rows are INFERRED (statistical test on burden summary stats).
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage1d")

URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Step 1: pull the disease list where any gene_burden evidence exists.
# Uses associatedDiseases with a datasources profile weighted to gene_burden so
# top rows are exactly those with a gene_burden contribution.
GQL_DISEASES = """
query D($ens: String!, $size: Int!) {
  target(ensemblId: $ens) {
    approvedSymbol
    associatedDiseases(
      page: {index: 0, size: $size}
      datasources: [{id: "gene_burden", weight: 1.0, propagate: true, required: true}]
    ) {
      count
      rows {
        disease { id name }
        datasourceScores { id score }
      }
    }
  }
}
"""

# Step 2: per (target, disease), pull the raw gene_burden evidence rows.
GQL_EV_ROWS = """
query EvRows($ens: [String!]!, $efo: String!) {
  disease(efoId: $efo) {
    id
    name
    evidences(ensemblIds: $ens, datasourceIds: ["gene_burden"], size: 100) {
      count
      rows {
        id
        datasourceId
        target { id approvedSymbol }
        disease { id name }
        pValueMantissa
        pValueExponent
        beta
        oddsRatio
        studyId
        studySampleSize
        studyCases
        ancestry
        ancestryId
        statisticalMethod
        statisticalMethodOverview
        literature
      }
    }
  }
}
"""


def _pvalue(row: dict) -> float | None:
    m = row.get("pValueMantissa")
    e = row.get("pValueExponent")
    if m is None or e is None:
        return None
    try:
        return float(m) * (10.0 ** float(e))
    except Exception:
        return None


def fetch_ot_gene_burden(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    stage1 = cfg["stage1"]
    if not stage1.get("associations", {}).get("ot_gene_burden", {}).get("enabled", False):
        log.info("ot_gene_burden disabled in config")
        return None

    p_max = float(stage1["associations"]["ot_gene_burden"].get("p_value_threshold", 1.0e-7))
    disease_size = int(stage1["associations"]["ot_gene_burden"].get("disease_scan_size", 100))
    sess = get_session()

    all_rows = []
    for sym, g in gene_ids.items():
        # Step 1: get diseases where gene_burden datasource contributed
        r = sess.post(URL, json={
            "query": GQL_DISEASES,
            "variables": {"ens": g.ensembl_id, "size": disease_size},
        }, timeout=60)
        if not r.ok:
            log.warning("OT associatedDiseases HTTP %d for %s | body=%s", r.status_code, sym, r.text[:400])
            continue
        js = r.json()
        if "errors" in js:
            log.warning("OT associatedDiseases GraphQL errors for %s | %s", sym, js["errors"][:2])
            continue
        t = (js.get("data") or {}).get("target")
        if not t:
            log.warning("OT: no target for %s", sym)
            continue

        rows = ((t.get("associatedDiseases") or {}).get("rows") or [])
        # Keep only rows where gene_burden actually contributed
        candidate_diseases = []
        for row in rows:
            dss = row.get("datasourceScores") or []
            gb = next((x for x in dss if x.get("id") == "gene_burden"), None)
            if gb and (gb.get("score") or 0) > 0:
                d = row.get("disease") or {}
                candidate_diseases.append((d.get("id"), d.get("name"), gb["score"]))
        log.info("OT gene_burden | %s | %d candidate diseases", sym, len(candidate_diseases))

        # Step 2: pull actual evidence rows per (target, disease)
        for efo, dname, score in candidate_diseases:
            r2 = sess.post(URL, json={
                "query": GQL_EV_ROWS,
                "variables": {"ens": [g.ensembl_id], "efo": efo},
            }, timeout=60)
            if not r2.ok:
                log.warning("OT evidences HTTP %d for %s/%s", r2.status_code, sym, efo)
                continue
            j2 = r2.json()
            if "errors" in j2:
                log.warning("OT evidences errors for %s/%s | %s", sym, efo, j2["errors"][:2])
                continue
            ev = ((j2.get("data") or {}).get("disease") or {}).get("evidences") or {}
            for row in (ev.get("rows") or []):
                p = _pvalue(row)
                if p is not None and p > p_max:
                    continue
                row["gene_symbol"] = sym
                row["ensembl_id"] = g.ensembl_id
                row["disease_id"] = efo
                row["disease_name"] = dname
                row["ot_gene_burden_ds_score"] = score
                row["pvalue"] = p
                row["certainty_tier"] = Tier.INFERRED.value
                row["source"] = SOURCES.get("ot_gene_burden", "opentargets_gene_burden_datasource")
                # Flatten
                row.pop("target", None)
                row.pop("disease", None)
                all_rows.append(row)

    if not all_rows:
        log.warning("OT gene_burden: no rows collected")
        empty = pd.DataFrame()
        safe_write_parquet(empty, outdir / "ot_gene_burden_evidences.parquet")
        return empty

    df = pd.DataFrame(all_rows)
    for col in ("literature", "allelicRequirements"):
        if col in df.columns:
            df[col] = df[col].apply(lambda x: ";".join(map(str, x)) if isinstance(x, list) else x)
    df = df.sort_values(["gene_symbol", "pvalue"]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "ot_gene_burden_evidences.parquet")
    log.info("ot_gene_burden_evidences.parquet written | rows=%d", len(df))
    return df
