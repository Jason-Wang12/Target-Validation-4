"""Fetch InterPro domain/family annotations for target UniProt IDs.

Cached to /mnt/results/rcvtc_runs/<run_id>/stage2/2c_uniprot/interpro_domains.parquet
"""
from pathlib import Path
import pandas as pd
import requests

UNIPROT_MAP = {
    "PCSK9": "Q8NBP7",
    "JAK2":  "O60674",
    "LDLR":  "P01130",
}

# Manual short-labels for prominent domains where InterPro names are verbose.
# Falls back to a truncated form of the InterPro name when not listed.
SHORT_LABEL = {
    "IPR000209": "Peptidase S8/S53",
    "IPR010259": "Inhibitor I9",
    "IPR041254": "CRD-1",
    "IPR041052": "CRD-2",
    "IPR041051": "CRD-3",
    "IPR000299": "FERM",
    "IPR000980": "SH2",
    "IPR000719": "Kinase (JH2/JH1)",
    "IPR001245": "Ser/Thr/Tyr kinase",
    "IPR000742": "EGF-like",
    "IPR001881": "EGF Ca-bind",
    "IPR049883": "EGF-like",
    "cd00112":   "LDLa",
    "cd00054":   "EGF_CA",
    "cd05078":   "JH2 pseudokinase",
    "cd10379":   "SH2",
    "cd13333":   "FERM C-lobe",
    "cd14205":   "JH1 kinase",
    "cd14473":   "FERM B-lobe",
    "cd04077":   "Peptidase S8",
    "cd16839":   "CRD",
    "IPR034193": "Proteinase K-like",
    "IPR050131": "Subtilisin-like",
    "IPR016251": "Jak/Tyk2 family",
    "IPR015500": "Subtilisin family",
    "IPR051221": "LDLR-related family",
}


def _shortlabel(acc: str, name: str) -> str:
    if acc in SHORT_LABEL:
        return SHORT_LABEL[acc]
    # Truncate to 20 chars
    return (name[:20] + "…") if len(name) > 20 else name


def fetch_interpro(gene: str, uniprot_id: str, timeout: int = 15) -> pd.DataFrame:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uniprot_id}/"
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
    r.raise_for_status()
    j = r.json()
    rows = []
    for entry in j.get("results", []):
        meta = entry.get("metadata", {})
        etype = meta.get("type", "").lower()
        if etype not in ("domain", "family"):
            continue
        acc = meta.get("accession", "")
        name = meta.get("name", "")
        source = meta.get("source_database", "")
        for prot in entry.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    start = frag.get("start")
                    end = frag.get("end")
                    if start is None or end is None:
                        continue
                    rows.append({
                        "gene_symbol": gene,
                        "uniprot_id": uniprot_id,
                        "interpro_id": acc,
                        "name": name,
                        "short_label": _shortlabel(acc, name),
                        "type": etype,
                        "source": source,
                        "start": int(start),
                        "end": int(end),
                    })
    return pd.DataFrame(rows)


def load_or_fetch(cache_path: Path, force: bool = False) -> pd.DataFrame:
    if cache_path.exists() and not force:
        return pd.read_parquet(cache_path)
    frames = []
    for gene, uni in UNIPROT_MAP.items():
        try:
            frames.append(fetch_interpro(gene, uni))
        except Exception as e:
            print(f"[interpro] FAILED {gene} {uni}: {e}")
    d = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    # Write via workspace-first (S3 no random-access), atomic mv
    tmp = Path("/workspace") / cache_path.name
    d.to_parquet(tmp, index=False)
    import shutil
    shutil.copy(tmp, cache_path)
    return d
