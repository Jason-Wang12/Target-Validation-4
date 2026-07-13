"""Stage 2d — AlphaFold per-residue pLDDT (ANNOTATED tier).

Approach:
    1. GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}
    2. Download the PDB file listed in the response
    3. Extract per-residue pLDDT from the B-factor column (AlphaFold convention)

Output: one row per (gene, residue) with the pLDDT value.
Downstream: variants can be joined by (gene, residue_number) to flag whether
they fall in structured (>=70) vs disordered (<50) regions.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from .certainty import Tier, SOURCES
from .utils.http_cache import get_session
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage2d")

AF_META = "https://alphafold.ebi.ac.uk/api/prediction/"


def _plddt_from_pdb(pdb_text: str) -> list[tuple[int, float]]:
    """Parse ATOM records; take one row per residue (the CA line)."""
    rows: dict[int, float] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        # PDB fixed-column layout: cols 22-26 residue number, cols 60-66 B-factor,
        # cols 13-16 atom name.
        try:
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            resi = int(line[22:26])
            bfac = float(line[60:66])
            rows[resi] = bfac
        except Exception:
            continue
    return sorted(rows.items())


def fetch_alphafold_plddt(cfg, gene_ids: dict, outdir: Path) -> Optional[pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    sess = get_session()

    all_rows = []
    for sym, g in gene_ids.items():
        if not g.uniprot_id:
            log.warning("skip %s: no UniProt ID", sym)
            continue
        r = sess.get(f"{AF_META}{g.uniprot_id}", timeout=30)
        if not r.ok:
            log.warning("AlphaFold meta failed | %s | %s", sym, r.status_code)
            continue
        d = r.json()
        if not isinstance(d, list) or not d:
            log.warning("AlphaFold no entry for %s", sym)
            continue
        entry = d[0]
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            log.warning("AlphaFold no pdbUrl for %s", sym)
            continue
        r2 = sess.get(pdb_url, timeout=60)
        if not r2.ok:
            log.warning("AlphaFold PDB download failed | %s | %s", sym, r2.status_code)
            continue
        residues = _plddt_from_pdb(r2.text)
        for resi, plddt in residues:
            all_rows.append({
                "gene_symbol": sym,
                "uniprot_id": g.uniprot_id,
                "residue_number": resi,
                "plddt": plddt,
                "certainty_tier": Tier.ANNOTATED.value,
                "source": SOURCES.get("alphafold", "alphafold_ebi_v6"),
            })
        log.info("  %s: %d residues (median pLDDT=%.1f)",
                 sym, len(residues),
                 sorted([p for _, p in residues])[len(residues) // 2] if residues else 0)

    if not all_rows:
        safe_write_parquet(pd.DataFrame(), outdir / "alphafold_plddt.parquet")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows).sort_values(["gene_symbol", "residue_number"]).reset_index(drop=True)
    safe_write_parquet(df, outdir / "alphafold_plddt.parquet")
    log.info("alphafold_plddt.parquet written | rows=%d", len(df))
    return df
