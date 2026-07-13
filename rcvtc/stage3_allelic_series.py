"""Stage 3 — Allelic series inference from Genebass mask cross-comparison.

Because Genebass provides gene-level burden summary stats (not individual-level
genotypes), the classical "functional-bin allelic series" test cannot be run
here. Instead, we use a mask-comparison proxy:

  For each (gene, phenotype) with p ≤ 1e-4 in Genebass, examine both masks:
    - pLoF (all high-confidence LoF variants collapsed)
    - missense|LC (missense + low-confidence LoF)

  Three canonical patterns emerge:

  LoF (haploinsufficiency): pLoF BETA has larger magnitude than missense|LC,
      and both are directionally concordant. Signature of dosage-sensitive genes
      where any null allele drives phenotype.

  GoF (activating): missense|LC BETA magnitude >= pLoF BETA magnitude and
      directions may diverge (LoF is protective, missense is deleterious).
      Signature of oncogenes and constitutively-active mutants.

  Recessive / LoF: pLoF BETA magnitude ≈ missense|LC BETA magnitude (both small),
      requires biallelic hits — burden on heterozygotes underpowered.

For each gene, we output all (phenotype, mask1_beta, mask2_beta, ratio,
concordance) rows for downstream mechanism verdict (stage 6).

Additionally, we annotate variants (from stage 1a gnomAD) with their
functional bin — using AlphaMissense (predicted tier) as the primary score
and MaveDB (measured, per-scoreset) where available. This bin field is
consumed by stage 9's evidence cards.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .certainty import Tier, SOURCES
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage3")


def _mask_comparison(gb: pd.DataFrame) -> pd.DataFrame:
    """Pivot Genebass into wide (mask side-by-side) and compute the proxy."""
    if gb.empty:
        return pd.DataFrame()

    # Genebass rows: one per (gene, annotation_mask, phenotype_slug)
    beta_col = "BETA_Burden" if "BETA_Burden" in gb.columns else "beta"
    se_col = "SE_Burden" if "SE_Burden" in gb.columns else "se"
    p_col = "Pvalue"

    keep_cols = ["gene", "phenotype_slug", "annotation_mask",
                 beta_col, se_col, p_col]
    keep_cols = [c for c in keep_cols if c in gb.columns]
    slim = gb[keep_cols].copy()

    # Bring in phenotype descriptor if present
    if "pheno_description" in gb.columns:
        slim["pheno_description"] = gb["pheno_description"]

    wide = slim.pivot_table(
        index=["gene", "phenotype_slug"],
        columns="annotation_mask",
        values=[beta_col, se_col, p_col],
        aggfunc="first",
    )
    # Flatten multi-index cols → "beta_pLoF", "se_pLoF", ...
    wide.columns = [f"{a}_{b}".replace("|", "_") for a, b in wide.columns]
    wide = wide.reset_index()

    # Compute derived fields where both masks are present
    lof = f"{beta_col}_pLoF"
    mis = f"{beta_col}_missense_LC"
    p_lof = f"{p_col}_pLoF"
    p_mis = f"{p_col}_missense_LC"
    if lof in wide and mis in wide:
        wide["pLoF_beta"] = wide[lof]
        wide["missense_beta"] = wide[mis]
        wide["ratio_missense_over_plof"] = wide[mis] / wide[lof]
        wide["direction_concordant"] = (
            (wide[lof] * wide[mis]) > 0
        )
        # Significance-weighted magnitude:  |β| × -log10(p)
        # Retained as auxiliary column even though the primary classifier is
        # β-magnitude based (see docstring).
        import numpy as np
        def _sig_mag(beta_col, pval_col):
            if beta_col in wide.columns and pval_col in wide.columns:
                p = wide[pval_col].astype(float).clip(lower=1e-320)
                return wide[beta_col].abs() * (-np.log10(p))
            elif beta_col in wide.columns:
                return wide[beta_col].abs()
            return pd.Series([np.nan]*len(wide), index=wide.index)
        wide["sig_mag_pLoF"] = _sig_mag(lof, p_lof)
        wide["sig_mag_missense"] = _sig_mag(mis, p_mis)

        # Pattern classification (β-magnitude based, with GoF-safe bins).
        #
        # Rationale: burden-test β is the per-carrier average effect within the
        # mask. LoF nulls typically have the largest per-carrier effect
        # (dosage), so |β_pLoF| >> |β_missense| is the classic haplo-
        # insufficiency fingerprint. When the missense β is comparable to or
        # exceeds pLoF, that pattern is unusual for pure LoF — it flags
        # dominant-negative / GoF / specific activating mechanisms.
        #
        # IMPORTANT LIMITATION: burden statistics can NOT distinguish a single
        # driver missense (e.g. JAK2 V617F) from a broad missense signal — the
        # per-carrier β gets diluted by non-driver missense carriers. Stage 6
        # supplements this with ClinVar single-variant / curated evidence.
        wide["mask_pattern"] = "unknown"
        both = wide[lof].notna() & wide[mis].notna()
        mag_lof = wide[lof].abs()
        mag_mis = wide[mis].abs()

        # LoF haploinsufficiency: pLoF magnitude clearly dominates (>=1.5×)
        cond_lof_haplo = both & (mag_lof >= mag_mis * 1.5) & wide["direction_concordant"]
        wide.loc[cond_lof_haplo, "mask_pattern"] = "lof_haploinsufficiency"

        # GoF / dominant-negative: missense >= pLoF magnitude
        cond_gof = both & (mag_mis > mag_lof) & wide["direction_concordant"]
        wide.loc[cond_gof, "mask_pattern"] = "gof_activating"

        # Comparable magnitudes with same direction (missense < pLoF < 1.5× missense)
        cond_mixed = both & wide["direction_concordant"] & (wide["mask_pattern"] == "unknown")
        wide.loc[cond_mixed, "mask_pattern"] = "mixed_mechanism"

        # Ambiguous: divergent directions
        cond_ambig = both & ~wide["direction_concordant"]
        wide.loc[cond_ambig, "mask_pattern"] = "direction_divergent"

        # Only one mask hit — can't discriminate
        wide.loc[wide[lof].isna() | wide[mis].isna(), "mask_pattern"] = "single_mask_only"

    return wide


def _annotate_variants_with_bins(cfg, outdir: Path) -> Optional[pd.DataFrame]:
    """Attach functional bin to each stage 1a variant (missense only)."""
    v_path = Path(cfg["outputs_dir"]) / "stage1" / "1a_variants" / "gnomad_variants.parquet"
    am_path = Path(cfg["outputs_dir"]) / "stage2" / "2b_alphamissense" / "alphamissense_scores.parquet"
    mv_path = Path(cfg["outputs_dir"]) / "stage2" / "2a_mavedb" / "mavedb_scores.parquet"

    if not v_path.exists():
        log.warning("Stage 1a variants missing — skipping bin annotation")
        return pd.DataFrame()

    variants = pd.read_parquet(v_path)
    miss = variants[variants["consequence"] == "missense_variant"].copy()

    # Attach AlphaMissense scores if available (predicted tier)
    if am_path.exists():
        am = pd.read_parquet(am_path)
        if not am.empty:
            miss = miss.merge(
                am[["gene_symbol", "variant_id", "am_class", "am_pathogenicity"]],
                on=["gene_symbol", "variant_id"], how="left",
            )
            log.info("bin-annotate | %d variants have AlphaMissense score",
                     miss["am_pathogenicity"].notna().sum())

    # Assign bin from am_pathogenicity (quartiles across all missense in the run)
    n_bins = cfg.get("stage3", {}).get("function_bins", {}).get("n_bins", 4)
    if "am_pathogenicity" in miss and miss["am_pathogenicity"].notna().any():
        # Rank across all genes together (comparable within pipeline)
        vals = miss.loc[miss["am_pathogenicity"].notna(), "am_pathogenicity"]
        bins = pd.qcut(vals, q=n_bins,
                       labels=[f"bin_{i+1}" for i in range(n_bins)])
        miss["am_bin"] = pd.Series(np.nan, index=miss.index, dtype="object")
        miss.loc[vals.index, "am_bin"] = bins.astype(str)
        miss["am_bin_source_tier"] = Tier.PREDICTED.value
    else:
        miss["am_bin"] = None
        miss["am_bin_source_tier"] = None

    # Overlay MaveDB (measured tier); when a variant has both, MaveDB overrides
    if mv_path.exists():
        mv = pd.read_parquet(mv_path)
        if not mv.empty and "hgvsp" in mv.columns:
            # MaveDB uses HGVS pro like "p.Gly2Leu"; gnomAD hgvsp is like "p.Gly2Leu"
            # Match by gene_symbol + hgvsp
            mv_slim = (mv.groupby(["gene_symbol", "hgvsp"])["score"]
                       .mean()
                       .reset_index()
                       .rename(columns={"score": "mavedb_score_mean"}))
            miss = miss.merge(mv_slim, on=["gene_symbol", "hgvsp"], how="left")
            log.info("bin-annotate | %d variants matched to MaveDB",
                     miss["mavedb_score_mean"].notna().sum())

            # MaveDB-based bins (per gene, quartiles) — override predicted where available
            for sym in miss["gene_symbol"].unique():
                sub_idx = miss.index[(miss["gene_symbol"] == sym) & miss["mavedb_score_mean"].notna()]
                if len(sub_idx) >= n_bins:
                    vals = miss.loc[sub_idx, "mavedb_score_mean"]
                    binlabs = pd.qcut(vals, q=n_bins, duplicates="drop",
                                      labels=False)
                    miss.loc[sub_idx, "mavedb_bin"] = [f"bin_{int(b)+1}" for b in binlabs]
                    miss.loc[sub_idx, "mavedb_bin_source_tier"] = Tier.MEASURED.value

    # Final "consensus_bin": prefer MaveDB (measured) over AlphaMissense (predicted)
    if "mavedb_bin" in miss:
        miss["consensus_bin"] = miss["mavedb_bin"].fillna(miss["am_bin"])
        miss["consensus_bin_tier"] = miss["mavedb_bin_source_tier"].fillna(miss["am_bin_source_tier"])
    else:
        miss["consensus_bin"] = miss["am_bin"]
        miss["consensus_bin_tier"] = miss["am_bin_source_tier"]

    safe_write_parquet(miss, outdir / "variants_binned.parquet")
    log.info("variants_binned.parquet written | rows=%d", len(miss))
    return miss


def run_stage3(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir = Path(cfg["outputs_dir"]) / "stage3"
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("stage3 start")

    counts = {}

    # 3a: Genebass mask comparison
    gb_path = Path(cfg["outputs_dir"]) / "stage1" / "1b_genebass" / "genebass_associations.parquet"
    if gb_path.exists():
        gb = pd.read_parquet(gb_path)
        wide = _mask_comparison(gb)
        safe_write_parquet(wide, outdir / "mask_comparison.parquet")
        log.info("mask_comparison.parquet written | rows=%d", len(wide))
        counts["mask_comparison_rows"] = len(wide)
        if not wide.empty and "mask_pattern" in wide.columns:
            log.info("pattern counts: %s", wide["mask_pattern"].value_counts().to_dict())
    else:
        log.warning("Genebass output missing — skipping mask comparison")
        counts["mask_comparison_rows"] = 0

    # 3b: Variant-level functional binning
    binned = _annotate_variants_with_bins(cfg, outdir)
    counts["variants_binned"] = len(binned) if binned is not None else 0

    log.info("stage3 counts: %s", counts)
    return {"counts": counts}
