"""
Certainty tier module — single source of truth for measured/inferred/predicted labels.

Every output field in the pipeline must be produced through one of these constructors.
Downstream synthesis (Stage 9) checks that no claim in an evidence card is missing a tier.

Tiers (per plan v3):
  MEASURED  — experimental data (MaveDB DMS scores, ClinVar PS3/BS3 functional evidence,
              gnomAD population allele frequencies as observed counts, GTEx median TPM,
              curated database entries from ClinVar/OMIM/HPO)
  INFERRED  — statistical test on summary statistics (Genebass burden, ExPheWas gene-level,
              OT gene_burden pass-through, FinnGen, tier-burden estimates)
  PREDICTED — computational model output (AlphaMissense, ESM, REVEL, LOEUF-anchored
              haploinsufficiency predictions, AlphaFold pLDDT)
  ANNOTATED — curator-verified structural / functional annotation (UniProt features,
              InterPro domains, PDB experimental interfaces). Treated as measured-adjacent
              for output purposes but tagged separately for provenance.

Usage:
    from rcvtc.certainty import Claim, Tier
    c = Claim(value=-0.039, tier=Tier.INFERRED, source="genebass_v450k_pLoF",
              note="annotation-mask aggregate", n=None)
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, field
from enum import Enum
from typing import Any, Optional


class Tier(str, Enum):
    MEASURED = "measured"
    INFERRED = "inferred"
    PREDICTED = "predicted"
    ANNOTATED = "annotated"


@dataclass
class Claim:
    """A single tier-tagged claim."""
    value: Any
    tier: Tier
    source: str                       # e.g. "genebass_v450k_pLoF", "alphamissense_v1"
    note: Optional[str] = None        # free-text caveat (e.g. "annotation-mask aggregate")
    n: Optional[int] = None           # sample / carrier size if applicable
    p_value: Optional[float] = None
    se: Optional[float] = None
    citation: Optional[str] = None    # DOI or PMID when the source is a specific paper

    def to_dict(self) -> dict:
        d = asdict(self)
        d["tier"] = self.tier.value  # serialize enum as string
        return d

    def to_display(self) -> str:
        """One-line display for markdown evidence cards."""
        v = self.value
        if isinstance(v, float):
            v = f"{v:.4g}"
        parts = [str(v)]
        detail = []
        if self.p_value is not None:
            detail.append(f"p={self.p_value:.2g}")
        if self.se is not None:
            detail.append(f"SE={self.se:.4g}")
        if self.n is not None:
            detail.append(f"n={self.n}")
        if detail:
            parts.append(f"({', '.join(detail)})")
        parts.append(f"[{self.tier.value} — {self.source}")
        if self.note:
            parts.append(f"; {self.note}")
        parts.append("]")
        return " ".join(parts).replace("[ ", "[").replace(" ; ", "; ")


@dataclass
class ClaimSet:
    """A group of related claims (e.g. all Genebass hits for one gene)."""
    claims: list[Claim] = field(default_factory=list)

    def add(self, claim: Claim) -> None:
        self.claims.append(claim)

    def to_list(self) -> list[dict]:
        return [c.to_dict() for c in self.claims]


def assert_tier_tagged(obj: Any, path: str = "") -> list[str]:
    """
    Walk a nested dict/list structure and return a list of dotted paths where a
    numeric or string value appears without an accompanying certainty tier.
    Used in Stage 9 to audit evidence-card outputs.

    Rules:
      - A dict with keys {'value', 'tier', 'source'} counts as a valid Claim shape.
      - A leaf value (int/float/str) is untagged unless its containing dict has a
        sibling '<key>_tier' or the whole dict is a Claim.
      - Lists are traversed by index.
    """
    unt = []
    if isinstance(obj, dict):
        if {"value", "tier", "source"}.issubset(obj.keys()):
            return unt
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                unt.extend(assert_tier_tagged(v, f"{path}.{k}" if path else k))
            elif isinstance(v, (int, float, str, bool)) and v is not None:
                # Metadata fields exempt from tagging
                if k in {"gene", "symbol", "ensembl_id", "uniprot_id", "hgnc_id",
                         "run_id", "phenotype", "annotation", "consequence",
                         "chrom", "pos", "ref", "alt", "hgvs_c", "hgvs_p",
                         "id", "name", "efo_id", "transcript_id", "canonical",
                         "phenotype_slug", "pheno_description", "outcome_id",
                         "outcome_label", "generated_at", "config_hash"}:
                    continue
                # Any *_tier field or *_source field is a companion, exempt
                if k.endswith("_tier") or k.endswith("_source") or k.endswith("_note"):
                    continue
                # Check sibling for its tier
                if f"{k}_tier" not in obj:
                    unt.append(f"{path}.{k}" if path else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list)):
                unt.extend(assert_tier_tagged(v, f"{path}[{i}]"))
    return unt


# Canonical source strings — used consistently across the pipeline
SOURCES = {
    "genebass": "genebass_ukb_wes_450k",
    "expheway": "expheway_v1",
    "ot_gene_burden": "opentargets_gene_burden_datasource",
    "finngen": "finngen_r12",
    "gnomad": "gnomad_v4_1",
    "alphamissense": "alphamissense_v1",
    "esm": "esm1v_t33_650M",
    "revel": "revel_dbnsfp_v4_5",
    "mavedb": "mavedb",
    "clinvar": "clinvar_2plus_star",
    "uniprot": "uniprot",
    "interpro": "interpro",
    "alphafold": "alphafold_db",
    "pdb": "rcsb_pdb",
    "gtex": "gtex_v8_median_tpm",
    "hpa": "hpa_v25",
    "opentargets": "opentargets_platform",
    "hpo": "hpo",
    "literature": "literature_search",
}
