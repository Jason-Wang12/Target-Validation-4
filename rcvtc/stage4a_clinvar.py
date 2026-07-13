"""Stage 4a — ClinVar pathogenic/benign variants (MEASURED-adjacent, ANNOTATED tier).

Uses NCBI E-utils:
    esearch clinvar term='LDLR[gene] AND (Pathogenic/Likely_pathogenic[CLIN])'
    esummary clinvar id=<batch of UIDs>

Note: the esearch [CLIN] filter is a broad match, so we re-filter on the
`germline_classification.description` field in the summary record. The
`review_status` field is mapped to ClinVar star tier for downstream filtering.

Configurable: cfg.stage4.clinvar.min_stars (default 2, "criteria provided,
multiple submitters no conflicts" or better).
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

log = logging.getLogger("rcvtc.stage4a")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Review-status → star mapping (public ClinVar rubric)
REVIEW_STAR = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,   # legacy label
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no interpretation for the single variant": 0,
}

# Which germline_classification.description strings are "pathogenic-ish"
PATH_LABELS = {
    "pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
    "pathogenic, low penetrance", "likely pathogenic, low penetrance",
}
BENIGN_LABELS = {
    "benign", "likely benign", "benign/likely benign",
}


def _eutils_get(sess, endpoint: str, params: dict, max_retries: int = 5) -> Optional[dict]:
    """GET with 429/5xx-aware retry and pacing (NCBI: 3 req/s anonymous)."""
    for attempt in range(max_retries):
        r = sess.get(f"{EUTILS}/{endpoint}", params=params, timeout=60)
        if r.ok:
            return r.json()
        if r.status_code in (429, 502, 503):
            # Honor Retry-After if present, else exponential backoff from 1s
            ra = r.headers.get("Retry-After")
            wait = int(ra) if ra and ra.isdigit() else min(2 ** attempt, 8)
            log.debug("eutils %s %s → sleep %ds", endpoint, r.status_code, wait)
            time.sleep(wait)
            continue
        log.warning("eutils %s failed | status=%d | body=%s",
                    endpoint, r.status_code, r.text[:200])
        return None
    log.warning("eutils %s: retries exhausted", endpoint)
    return None


def _fetch_variant_ids(sess, sym: str, sig_label: str) -> list[str]:
    """Return list of ClinVar UIDs for the given gene/significance combo."""
    j = _eutils_get(sess, "esearch.fcgi", {
        "db": "clinvar",
        "term": f"({sym}[gene]) AND ({sig_label}[Clinical_significance])",
        "retmax": 5000,
        "retmode": "json",
    })
    time.sleep(0.34)   # 3 req/s pacing
    if not j:
        return []
    return j.get("esearchresult", {}).get("idlist", [])


def _fetch_summaries(sess, uids: list[str], batch_size: int = 100) -> list[dict]:
    """Fetch esummary for a batch of UIDs (batches of 100 = NCBI recommended)."""
    out = []
    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i + batch_size]
        j = _eutils_get(sess, "esummary.fcgi", {
            "db": "clinvar", "id": ",".join(chunk), "retmode": "json"
        })
        if not j:
            continue
        for uid in j.get("result", {}).get("uids", []):
            rec = j["result"].get(uid) or {}
            rec["uid"] = uid
            out.append(rec)
        time.sleep(0.34)  # 3 req/s pacing
    return out


def _rec_to_row(rec: dict, sym: str) -> Optional[dict]:
    """Flatten one ClinVar esummary record into a row; None if unusable."""
    gcls = rec.get("germline_classification") or {}
    label = (gcls.get("description") or "").strip().lower()
    if not label:
        return None
    review = (gcls.get("review_status") or "").strip().lower()
    stars = REVIEW_STAR.get(review, 0)

    if label in PATH_LABELS:
        clinsig_bucket = "pathogenic_lp"
    elif label in BENIGN_LABELS:
        clinsig_bucket = "benign_lb"
    else:
        clinsig_bucket = "other"

    vs_list = rec.get("variation_set") or []
    vs = vs_list[0] if vs_list else {}
    locs = vs.get("variation_loc") or []
    grch38 = next((l for l in locs if l.get("assembly_name") == "GRCh38"), None)
    chrom = grch38.get("chr") if grch38 else None
    pos = grch38.get("start") if grch38 else None
    variant_name = vs.get("variation_name") or rec.get("title")

    traits = [
        t.get("trait_name")
        for t in (gcls.get("trait_set") or [])
        if t.get("trait_name") and t.get("trait_name") != "not provided"
    ]

    return {
        "gene_symbol": sym,
        "clinvar_uid": rec.get("uid"),
        "accession": rec.get("accession"),
        "variant_name": variant_name,
        "cdna_change": vs.get("cdna_change"),
        "protein_change": rec.get("protein_change"),
        "variant_type": vs.get("variant_type"),
        "chrom": chrom,
        "pos_grch38": pos,
        "canonical_spdi": vs.get("canonical_spdi"),
        "germline_classification": gcls.get("description"),
        "review_status": gcls.get("review_status"),
        "stars": stars,
        "clinsig_bucket": clinsig_bucket,
        "last_evaluated": gcls.get("last_evaluated"),
        "traits": ";".join(traits) if traits else None,
        "certainty_tier": Tier.ANNOTATED.value,   # curated per-variant classification
        "source": SOURCES.get("clinvar", "clinvar"),
    }


def fetch_clinvar_variants(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()
    min_stars = cfg.get("stage4", {}).get("clinvar", {}).get("min_stars", 2)

    rows = []
    for sym in gene_ids.keys():
        # Two searches: pathogenic side and benign side (kept separate for logging)
        for sig_label in ("Pathogenic/Likely_pathogenic", "Benign/Likely_benign"):
            uids = _fetch_variant_ids(sess, sym, sig_label)
            if not uids:
                log.info("ClinVar | %s | %s | 0 hits", sym, sig_label)
                continue
            log.info("ClinVar | %s | %s | %d uids", sym, sig_label, len(uids))
            recs = _fetch_summaries(sess, uids)
            gene_rows = [r for rec in recs if (r := _rec_to_row(rec, sym))]
            rows.extend(gene_rows)
            log.info("  %s | %s | %d records flattened", sym, sig_label, len(gene_rows))

    if not rows:
        empty = pd.DataFrame()
        safe_write_parquet(empty, outdir / "clinvar_variants.parquet")
        return empty

    df = pd.DataFrame(rows).drop_duplicates(subset=["clinvar_uid"]).reset_index(drop=True)
    # Report before filtering
    log.info(
        "ClinVar | %d unique variants across %d genes | stars>=%d passes %d/%d",
        len(df), df["gene_symbol"].nunique(), min_stars,
        (df["stars"] >= min_stars).sum(), len(df),
    )

    # Keep only records at or above the star threshold + real bucket
    keep = df[(df["stars"] >= min_stars) & (df["clinsig_bucket"] != "other")].copy()
    safe_write_parquet(keep, outdir / "clinvar_variants.parquet")
    log.info("clinvar_variants.parquet written | rows=%d", len(keep))
    return keep
