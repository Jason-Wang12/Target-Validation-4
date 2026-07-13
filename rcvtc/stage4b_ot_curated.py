"""Stage 4b — Open Targets curated evidence (ANNOTATED tier).

Retrieves gene-disease evidence rows from the "curated" datasources:
    eva                    — ClinVar-derived, curated by EBI
    orphanet               — rare-disease curation
    genomics_england       — GEL PanelApp
    gene2phenotype          — DECIPHER
    uniprot_literature     — UniProt curated disease associations

Uses the same two-step flow as stage1d (associatedDiseases → per-disease evidences)
but keys on `datasourceIds` list rather than a single datasource.
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

log = logging.getLogger("rcvtc.stage4b")

OT_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Step 1: get diseases associated by any of the curated sources
QUERY_DISEASES = """
query D($ens: String!, $srcs: [DatasourceSettingsInput!]!, $size: Int!) {
  target(ensemblId: $ens) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: $size}, datasources: $srcs) {
      count
      rows {
        disease { id name therapeuticAreas { id name } }
        score
        datasourceScores { id score }
      }
    }
  }
}
"""

# Step 2: per (target, disease), pull curated evidence
QUERY_EVIDENCES = """
query Ev($ens: [String!]!, $efo: String!, $srcIds: [String!]!, $size: Int!) {
  disease(efoId: $efo) {
    evidences(ensemblIds: $ens, datasourceIds: $srcIds, size: $size) {
      count
      rows {
        id
        datasourceId
        datatypeId
        score
        variantRsId
        variantAminoacidDescriptions
        variantFunctionalConsequence { id label }
        confidence
        clinicalSignificances
        allelicRequirements
        cohortId
        studyId
        literature
        diseaseFromSource
        diseaseFromSourceId
      }
    }
  }
}
"""


def fetch_ot_curated(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()

    src_ids = cfg["stage4"]["ot_curated"]["datasource_ids"]
    max_ev = cfg["stage4"]["ot_curated"].get("max_evidences_per_disease", 200)

    # associatedDiseases uses the DatasourceSettingsInput type (same as stage1d)
    src_settings = [{"id": s, "weight": 1.0, "propagate": True, "required": True}
                    for s in src_ids]

    diseases_rows = []
    ev_rows = []

    for sym, g in gene_ids.items():
        # Step 1: associated diseases
        r = sess.post(
            OT_GQL,
            json={"query": QUERY_DISEASES,
                  "variables": {"ens": g.ensembl_id,
                                "srcs": src_settings,
                                "size": 500}},
            timeout=60,
        )
        if not r.ok:
            log.warning("OT curated D query failed | %s | %s | body=%s",
                        sym, r.status_code, r.text[:400])
            continue
        j = r.json()
        target = ((j.get("data") or {}).get("target") or {})
        assoc = (target.get("associatedDiseases") or {})
        rows = assoc.get("rows") or []
        log.info("OT curated | %s | %d diseases (of count=%s)",
                 sym, len(rows), assoc.get("count"))

        for row in rows:
            d = row.get("disease", {})
            diseases_rows.append({
                "gene_symbol": sym,
                "ensembl_id": g.ensembl_id,
                "disease_id": d.get("id"),
                "disease_name": d.get("name"),
                "therapeutic_areas": ";".join(
                    ta.get("name") for ta in (d.get("therapeuticAreas") or [])
                ) or None,
                "overall_score": row.get("score"),
                "datasource_scores": ";".join(
                    f"{ds.get('id')}={ds.get('score'):.3f}"
                    for ds in (row.get("datasourceScores") or []) if ds.get("id") in src_ids
                ),
                "certainty_tier": Tier.ANNOTATED.value,
                "source": SOURCES.get("opentargets", "opentargets_platform"),
            })

        # Step 2: per-disease evidence pull (only for the top-N of the shortlist to keep runtime bounded)
        for row in rows:
            efo = row.get("disease", {}).get("id")
            if not efo:
                continue
            r2 = sess.post(
                OT_GQL,
                json={"query": QUERY_EVIDENCES,
                      "variables": {"ens": [g.ensembl_id],
                                    "efo": efo,
                                    "srcIds": src_ids,
                                    "size": max_ev}},
                timeout=60,
            )
            if not r2.ok:
                log.warning("OT curated Ev query failed | %s | %s | %s",
                            sym, efo, r2.status_code)
                continue
            j2 = r2.json()
            disease_obj = ((j2.get("data") or {}).get("disease") or {})
            ev = (disease_obj.get("evidences") or {}).get("rows") or []
            for e in ev:
                vfc = (e.get("variantFunctionalConsequence") or {})
                ev_rows.append({
                    "gene_symbol": sym,
                    "ensembl_id": g.ensembl_id,
                    "disease_id": efo,
                    "disease_name": row.get("disease", {}).get("name"),
                    "evidence_id": e.get("id"),
                    "datasource_id": e.get("datasourceId"),
                    "datatype_id": e.get("datatypeId"),
                    "score": e.get("score"),
                    "variant_rs_id": e.get("variantRsId"),
                    "variant_aa_descriptions": ";".join(e.get("variantAminoacidDescriptions") or []) or None,
                    "variant_functional_consequence_id": vfc.get("id"),
                    "variant_functional_consequence_label": vfc.get("label"),
                    "confidence": e.get("confidence"),
                    "clinical_significances": ";".join(e.get("clinicalSignificances") or []) or None,
                    "allelic_requirements": ";".join(e.get("allelicRequirements") or []) or None,
                    "cohort_id": e.get("cohortId"),
                    "study_id": e.get("studyId"),
                    "literature": ";".join(e.get("literature") or []) or None,
                    "disease_from_source": e.get("diseaseFromSource"),
                    "certainty_tier": Tier.ANNOTATED.value,
                    "source": SOURCES.get("opentargets", "opentargets_platform"),
                })
            time.sleep(0.05)   # gentle pacing across OT
        log.info("  %s | %d curated evidence rows so far", sym, sum(1 for e in ev_rows if e['gene_symbol'] == sym))

    diseases_df = pd.DataFrame(diseases_rows)
    ev_df = pd.DataFrame(ev_rows)

    safe_write_parquet(diseases_df, outdir / "ot_curated_diseases.parquet")
    safe_write_parquet(ev_df, outdir / "ot_curated_evidences.parquet")
    log.info("ot_curated: %d diseases, %d evidences written", len(diseases_df), len(ev_df))
    return ev_df
