"""HGNC ↔ Ensembl ↔ UniProt resolution via Open Targets + UniProt."""
from __future__ import annotations
import logging
from typing import Optional
from dataclasses import dataclass, asdict

from .http_cache import get_session

log = logging.getLogger("rcvtc.ids")

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"


@dataclass
class GeneIDs:
    symbol: str                  # HGNC-approved symbol
    ensembl_id: str              # ENSG...
    uniprot_id: Optional[str]    # canonical Swiss-Prot accession
    canonical_transcript_id: Optional[str] = None   # ENST... (Ensembl canonical)

    def to_dict(self) -> dict:
        return asdict(self)


def resolve_gene(symbol_or_id: str) -> Optional[GeneIDs]:
    """Resolve a symbol or Ensembl ID to a GeneIDs record via Open Targets."""
    sess = get_session()
    query = symbol_or_id.strip()
    is_ensg = query.upper().startswith("ENSG")

    if is_ensg:
        # Direct target query
        gql = """
        query T($id: String!) {
          target(ensemblId: $id) {
            id approvedSymbol
            proteinIds { id source }
            canonicalTranscript { id }
          }
        }"""
        variables = {"id": query.upper()}
        r = sess.post(OT_URL, json={"query": gql, "variables": variables}, timeout=30)
    else:
        # Search-then-fetch
        gql = """
        query S($q: String!) {
          search(queryString: $q, entityNames: ["target"]) {
            hits { id name entity }
          }
        }"""
        r = sess.post(OT_URL, json={"query": gql, "variables": {"q": query}}, timeout=30)

    if not r.ok:
        log.warning("Open Targets returned HTTP %s for %s", r.status_code, symbol_or_id)
        return None
    data = r.json().get("data") or {}

    if is_ensg:
        t = data.get("target")
        if not t:
            log.warning("Open Targets returned no target for %s", symbol_or_id)
            return None
        return _pack(t)

    hits = (data.get("search") or {}).get("hits") or []
    hits = [h for h in hits if h.get("entity") == "target"]
    if not hits:
        log.warning("Open Targets search returned no target hits for %s", symbol_or_id)
        return None

    # Prefer exact symbol match; else first
    best = None
    q_upper = query.upper()
    for h in hits:
        if h["name"].upper() == q_upper:
            best = h
            break
    if best is None:
        best = hits[0]

    # Follow-up fetch to get UniProt + canonical transcript
    gql2 = """
    query T($id: String!) {
      target(ensemblId: $id) {
        id approvedSymbol
        proteinIds { id source }
        canonicalTranscript { id }
      }
    }"""
    r2 = sess.post(OT_URL, json={"query": gql2, "variables": {"id": best["id"]}}, timeout=30)
    if not r2.ok:
        return GeneIDs(symbol=best["name"], ensembl_id=best["id"], uniprot_id=None)
    t = (r2.json().get("data") or {}).get("target")
    if not t:
        return GeneIDs(symbol=best["name"], ensembl_id=best["id"], uniprot_id=None)
    return _pack(t)


def _pack(t: dict) -> GeneIDs:
    ens = t["id"]
    sym = t["approvedSymbol"]
    uniprot = None
    for pid in t.get("proteinIds") or []:
        # Prefer uniprot_swissprot; else any uniprot
        if pid.get("source") == "uniprot_swissprot":
            uniprot = pid["id"]
            break
    if uniprot is None:
        for pid in t.get("proteinIds") or []:
            if "uniprot" in (pid.get("source") or "").lower():
                uniprot = pid["id"]
                break
    canon = None
    ct = t.get("canonicalTranscript")
    if ct:
        canon = ct.get("id")
    return GeneIDs(symbol=sym, ensembl_id=ens, uniprot_id=uniprot, canonical_transcript_id=canon)


def resolve_many(symbols: list[str]) -> dict[str, Optional[GeneIDs]]:
    """Resolve a batch of gene symbols to IDs. Preserves input order in the result dict."""
    out: dict[str, Optional[GeneIDs]] = {}
    for s in symbols:
        out[s] = resolve_gene(s)
        if out[s] is None:
            log.warning("Could not resolve gene: %s", s)
    return out
