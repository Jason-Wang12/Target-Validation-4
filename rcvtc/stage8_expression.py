"""Stage 8 — Expression / tissue context (Human Protein Atlas consensus tissue nTPM).

Fetches per-tissue nTPM from HPA per gene (via https://www.proteinatlas.org/{Ensembl}.xml).
Uses the HPA `consensusTissue` assay (51 tissues; harmonized HPA+GTEx+FANTOM).
Also emits tissue-enrichment classification and a top-N ranking.

All rows are MEASURED tier — direct RNA-seq nTPM values from HPA v25.1.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd
from lxml import etree

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage8")

HPA_XML_TMPL = "https://www.proteinatlas.org/{ensg}.xml"


def _parse_hpa_xml(xml_bytes: bytes) -> dict:
    """Extract consensusTissue nTPM + tissue-specificity summary from an HPA gene XML."""
    tree = etree.fromstring(xml_bytes)

    tissues: dict[str, dict] = {}
    specificity = None
    tissue_enriched_in = None
    specificity_score = None
    for rna in tree.findall(".//rnaExpression"):
        if rna.get("assayType") != "consensusTissue":
            continue
        for d in rna.findall("data"):
            t_el = d.find("tissue")
            if t_el is None or t_el.text is None:
                continue
            tissue = t_el.text
            organ = t_el.get("organ")
            ntpm_el = d.find("level[@type='normalizedRNAExpression']")
            ntpm = None
            if ntpm_el is not None:
                try:
                    ntpm = float(ntpm_el.get("expRNA"))
                except (TypeError, ValueError):
                    ntpm = None
            tissues[tissue] = {"organ": organ, "ntpm": ntpm}
        # tissue-specificity summary for the same assay
        rs = rna.find("rnaSpecificity")
        if rs is not None:
            specificity = rs.get("specificity")
            enriched = [t.text for t in rs.findall("tissue") if t.text]
            if enriched:
                tissue_enriched_in = ";".join(enriched)
        # specificity score (integer 0-10)
        ss = rna.find("specificityScore")
        if ss is not None and ss.text:
            try:
                specificity_score = float(ss.text)
            except (TypeError, ValueError):
                pass
        break  # consensusTissue is unique; stop

    return {
        "tissues": tissues,
        "specificity": specificity,
        "tissue_enriched_in": tissue_enriched_in,
        "specificity_score": specificity_score,
    }


def run_stage8(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage8"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    sess = get_session()

    log.info("stage8 start | %d genes", len(gene_ids))

    long_rows = []
    summary_rows = []
    for sym, ids in gene_ids.items():
        # ids is a GeneIDs dataclass (or None if unresolved)
        ensg = getattr(ids, "ensembl_id", None) if ids is not None else None
        if not ensg:
            log.warning("no Ensembl id for %s", sym)
            continue
        url = HPA_XML_TMPL.format(ensg=ensg)
        r = sess.get(url, timeout=45)
        if not r.ok:
            log.warning("HPA XML fetch failed | %s | %s | %s", sym, ensg, r.status_code)
            continue
        try:
            data = _parse_hpa_xml(r.content)
        except Exception as e:
            log.warning("HPA XML parse failed | %s | %s", sym, e)
            continue

        # Long-format rows
        for tissue, tinfo in data["tissues"].items():
            long_rows.append({
                "gene_symbol": sym,
                "ensembl_id": ensg,
                "tissue": tissue,
                "organ": tinfo.get("organ"),
                "ntpm": tinfo.get("ntpm"),
                "certainty_tier": Tier.MEASURED.value,
                "source": SOURCES.get("hpa", "hpa_v25"),
            })

        # Summary + top tissues
        ntpm_series = pd.Series({t: v.get("ntpm") for t, v in data["tissues"].items()})
        ntpm_series = ntpm_series.dropna()
        top_tissues = ntpm_series.sort_values(ascending=False).head(5)
        summary_rows.append({
            "gene_symbol": sym,
            "ensembl_id": ensg,
            "n_tissues": len(ntpm_series),
            "max_ntpm": float(ntpm_series.max()) if len(ntpm_series) else None,
            "max_tissue": ntpm_series.idxmax() if len(ntpm_series) else None,
            "median_ntpm": float(ntpm_series.median()) if len(ntpm_series) else None,
            "tissue_specificity": data["specificity"],
            "tissue_enriched_in": data["tissue_enriched_in"],
            "specificity_score": data["specificity_score"],
            "top5_tissues": ";".join(f"{t}:{v:.1f}" for t, v in top_tissues.items()),
            "certainty_tier": Tier.MEASURED.value,
            "source": SOURCES.get("hpa", "hpa_v25"),
        })
        log.info(
            "hpa | %s | top=%s(%.1f) | specificity=%s",
            sym,
            (ntpm_series.idxmax() if len(ntpm_series) else "-"),
            (ntpm_series.max() if len(ntpm_series) else float("nan")),
            data["specificity"],
        )

    long_df = pd.DataFrame(long_rows)
    sum_df = pd.DataFrame(summary_rows)
    safe_write_parquet(long_df, outdir / "hpa_tissue_ntpm.parquet")
    safe_write_parquet(sum_df, outdir / "hpa_expression_summary.parquet")
    log.info("stage8 wrote | %d long-rows | %d summary-rows", len(long_df), len(sum_df))
    return {"counts": {"hpa_tissue_ntpm": len(long_df), "hpa_summary": len(sum_df)}}
