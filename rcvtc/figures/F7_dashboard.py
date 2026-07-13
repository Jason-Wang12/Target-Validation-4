"""F7 — Target validation dashboard composite.

Landscape 16:9 slide. Rows = 3 genes (PCSK9, JAK2, LDLR).
Columns:
  1) Header card (gene name, verdict badge, LOEUF, rank).
  2) Genebass burden — dual-mask βs and p-values (mini table).
  3) ClinVar 2⁺★ pathogenic — count + missense/LoF/other breakdown, top variants.
  4) AlphaMissense summary — n likely_pathogenic + mean pathogenicity.
  5) HPA top-tissue enrichment — top-3 tissues with nTPM bars.
  6) Top curated disease — top-3 rows with weighted_sum.

All values are pulled from source parquets — no hard-coded numbers.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from figures.palette import GENE_COLOR, VERDICT_COLOR, AM_CLASS_COLOR
from figures.style import apply_style, save

GENES = ["PCSK9", "JAK2", "LDLR"]

_LOF_TOKENS = ("*", "ter", "fs", "del", "dup")


def _classify_protein_change(pc) -> str:
    if pc is None or (isinstance(pc, float) and np.isnan(pc)):
        return "unknown"
    s = str(pc).lower()
    tok = s.split(",")[0].strip()
    if any(t in tok for t in _LOF_TOKENS):
        return "lof"
    if re.match(r"[a-z]+\d+[a-z*]", tok):
        return "missense"
    return "unknown"


def _fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and (np.isnan(p) or p <= 0)):
        return "n/a"
    if p < 1e-100:
        return "<1e-100"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3g}"


def _pretty_pheno(slug: str) -> str:
    """Turn Genebass slug into readable label."""
    mapping = {
        "ldl-direct": "LDL-C direct",
        "cholesterol": "Cholesterol",
        "apolipoprotein-b": "Apolipoprotein B",
        "platelet-count": "Platelet count",
        "platelet-crit": "Platelet crit",
        "eosinophill-count": "Eosinophil count",
        "igf-1": "IGF-1",
        "date-e78-first-reported-disorders-of-lipoprotein-metabolism-and-other-lipidaemias": "Lipid disorder (E78)",
    }
    return mapping.get(slug, slug.replace("-", " ").capitalize())


def _draw_header(ax, gene: str, verdict: str, rank: int, loeuf: float, composite: float):
    ax.axis("off")
    color = GENE_COLOR[gene]
    verdict_c = VERDICT_COLOR.get(verdict, "#666666")
    ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.15,
                                     transform=ax.transAxes))
    ax.text(0.05, 0.82, gene, fontsize=22, fontweight="bold",
            color=color, transform=ax.transAxes)
    ax.text(0.05, 0.68, f"Rank #{rank}", fontsize=11, color="#555",
            transform=ax.transAxes)
    # Verdict pill
    label_map = {"lof_haploinsufficiency": "LoF · haploinsufficiency",
                 "gof_activating": "GoF · activating",
                 "unknown": "unresolved",
                 "mixed_mechanism": "mixed mechanism"}
    verdict_label = label_map.get(verdict, verdict)
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.05, 0.42), 0.85, 0.16,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        transform=ax.transAxes,
        facecolor=verdict_c, alpha=0.85, edgecolor="none",
    ))
    ax.text(0.475, 0.51, verdict_label, fontsize=10.5, fontweight="bold",
            color="white", ha="center", va="center", transform=ax.transAxes)
    # LOEUF & composite
    ax.text(0.05, 0.28, f"LOEUF: {loeuf:.2f}", fontsize=9, transform=ax.transAxes)
    ax.text(0.05, 0.18, f"Composite score: {composite:.1f}", fontsize=9, transform=ax.transAxes)


def _draw_genebass_card(ax, gb: pd.DataFrame, mc: pd.DataFrame, gene: str):
    ax.axis("off")
    ax.text(0.02, 0.95, "Genebass burden (top phenotype)",
            fontsize=10, fontweight="semibold", transform=ax.transAxes, va="top")

    # Find top phenotype: pick one with strongest signal in either mask
    gb_g = gb[gb["gene"] == gene]
    if len(gb_g) == 0:
        ax.text(0.02, 0.7, "no signal", fontsize=9, color="#888", transform=ax.transAxes)
        return
    idxmin = gb_g["Pvalue_Burden"].idxmin()
    top = gb_g.loc[idxmin]
    pheno_slug = top["phenotype_slug"]
    label = _pretty_pheno(pheno_slug)

    ax.text(0.02, 0.80, f"Phenotype: {label}",
            fontsize=9, transform=ax.transAxes, va="top", color="#333")

    # Pull both masks
    both = mc[(mc["gene"] == gene) & (mc["phenotype_slug"] == pheno_slug)]
    if len(both):
        r = both.iloc[0]
        beta_lof = r["BETA_Burden_pLoF"]
        beta_mis = r["BETA_Burden_missense_LC"]
        p_lof = r["Pvalue_pLoF"]
        p_mis = r["Pvalue_missense_LC"]
    else:
        # Fallback: use whichever mask this row was
        beta_lof = top["BETA_Burden"] if top["annotation_mask"] == "pLoF" else None
        beta_mis = top["BETA_Burden"] if top["annotation_mask"] == "missense|LC" else None
        p_lof = top["Pvalue_Burden"] if top["annotation_mask"] == "pLoF" else None
        p_mis = top["Pvalue_Burden"] if top["annotation_mask"] == "missense|LC" else None

    def _row(y, mask_label, beta, p):
        ax.text(0.02, y, mask_label, fontsize=8, color="#555",
                transform=ax.transAxes, va="top")
        if beta is None or (isinstance(beta, float) and np.isnan(beta)):
            ax.text(0.30, y, "n/a", fontsize=8.5, transform=ax.transAxes, va="top", color="#888")
        else:
            ax.text(0.30, y, f"β = {beta:+.3f}",
                    fontsize=8.5, transform=ax.transAxes, va="top", color="#111")
            ax.text(0.60, y, f"p = {_fmt_p(p)}",
                    fontsize=8.5, transform=ax.transAxes, va="top", color="#111")

    _row(0.60, "pLoF mask", beta_lof, p_lof)
    _row(0.44, "missense|LC", beta_mis, p_mis)


def _draw_clinvar_card(ax, cv: pd.DataFrame, gene: str):
    ax.axis("off")
    ax.text(0.02, 0.95, "ClinVar 2⁺★ pathogenic",
            fontsize=10, fontweight="semibold", transform=ax.transAxes, va="top")
    g = cv[(cv["gene_symbol"] == gene) & (cv["clinsig_bucket"] == "pathogenic_lp")]
    n_total = len(g)
    n_miss = 0
    n_lof = 0
    for _, r in g.iterrows():
        cls = _classify_protein_change(r["protein_change"])
        if cls == "missense":
            n_miss += 1
        elif cls == "lof":
            n_lof += 1
    ax.text(0.02, 0.80, f"n = {n_total}", fontsize=11, fontweight="bold",
            transform=ax.transAxes, va="top")
    ax.text(0.02, 0.66, f"missense: {n_miss}   LoF: {n_lof}   other: {n_total - n_miss - n_lof}",
            fontsize=8.5, transform=ax.transAxes, va="top", color="#333")

    # Top variants: for high-count LDLR, show one representative missense + one LoF token
    # For low-count PCSK9/JAK2, list up to 3 protein_change tokens
    def _first_token(pc):
        if pc is None or (isinstance(pc, float) and np.isnan(pc)):
            return None
        return str(pc).split(",")[0].strip()

    tokens = [(_first_token(r["protein_change"]), _classify_protein_change(r["protein_change"]))
              for _, r in g.iterrows()]
    tokens = [(t, c) for (t, c) in tokens if t]
    if len(tokens) <= 5:
        ax.text(0.02, 0.50, "Variants:", fontsize=8, color="#555",
                transform=ax.transAxes, va="top")
        y0 = 0.38
        for tok, cls in tokens[:5]:
            ax.text(0.02, y0, tok, fontsize=8, color="#111", transform=ax.transAxes, va="top")
            y0 -= 0.09
    else:
        # LDLR case — enumerate hotspot residues (top-5 most common protein_change tokens)
        cnt = pd.Series([t for t, _ in tokens]).value_counts().head(5)
        ax.text(0.02, 0.50, "Hotspot variants:", fontsize=8, color="#555",
                transform=ax.transAxes, va="top")
        y0 = 0.38
        for tok, c in cnt.items():
            ax.text(0.02, y0, f"{tok}  (×{c})", fontsize=8,
                    color="#111", transform=ax.transAxes, va="top")
            y0 -= 0.09


def _draw_am_card(ax, am: pd.DataFrame, gene: str):
    ax.axis("off")
    ax.text(0.02, 0.95, "AlphaMissense",
            fontsize=10, fontweight="semibold", transform=ax.transAxes, va="top")
    g = am[am["gene_symbol"] == gene]
    n = len(g)
    n_path = int((g["am_class"] == "likely_pathogenic").sum())
    n_amb = int((g["am_class"] == "ambiguous").sum())
    n_ben = int((g["am_class"] == "likely_benign").sum())
    mean_p = float(g["am_pathogenicity"].mean()) if n else float("nan")

    ax.text(0.02, 0.80, f"n scored: {n:,}", fontsize=9, transform=ax.transAxes, va="top")
    ax.text(0.02, 0.68, f"mean pathogenicity: {mean_p:.3f}",
            fontsize=9, transform=ax.transAxes, va="top")

    # Compact stacked bar (horizontal) showing class fractions
    if n:
        f_ben, f_amb, f_path = n_ben / n, n_amb / n, n_path / n
        x0, width = 0.02, 0.86
        y0, height = 0.36, 0.10
        for f, color, label in zip(
            [f_ben, f_amb, f_path],
            [AM_CLASS_COLOR["likely_benign"], AM_CLASS_COLOR["ambiguous"], AM_CLASS_COLOR["likely_pathogenic"]],
            ["benign", "amb.", "path."],
        ):
            ax.add_patch(mpatches.Rectangle((x0, y0), width * f, height,
                                             facecolor=color, edgecolor="none",
                                             transform=ax.transAxes))
            x0 += width * f
        # Counts below
        ax.text(0.02, 0.22,
                f"benign: {n_ben:,}   amb.: {n_amb:,}   likely_path.: {n_path:,}",
                fontsize=8, transform=ax.transAxes, va="top", color="#333")


def _draw_hpa_card(ax, hpa: pd.DataFrame, summary: pd.DataFrame, gene: str):
    ax.axis("off")
    s = summary[summary["gene_symbol"] == gene].iloc[0]
    ax.text(0.02, 0.95, "HPA tissue expression",
            fontsize=10, fontweight="semibold", transform=ax.transAxes, va="top")
    ax.text(0.02, 0.82, s["tissue_specificity"],
            fontsize=9, color="#333", transform=ax.transAxes, va="top")

    g = hpa[hpa["gene_symbol"] == gene].nlargest(3, "ntpm")
    y_start = 0.65
    max_ntpm = g["ntpm"].max() if len(g) else 1.0
    for i, (_, r) in enumerate(g.iterrows()):
        y = y_start - i * 0.20
        # tissue label
        ax.text(0.02, y + 0.02, r["tissue"], fontsize=8, color="#333",
                transform=ax.transAxes, va="top")
        # bar
        bar_len = 0.5 * (r["ntpm"] / max_ntpm)
        ax.add_patch(mpatches.Rectangle((0.40, y - 0.05), bar_len, 0.06,
                                         facecolor=GENE_COLOR[gene], alpha=0.85,
                                         transform=ax.transAxes, edgecolor="none"))
        ax.text(0.40 + bar_len + 0.01, y - 0.02, f"{r['ntpm']:.1f}",
                fontsize=8, color="#111", transform=ax.transAxes, va="top")


def _draw_disease_card(ax, dz: pd.DataFrame, gene: str):
    ax.axis("off")
    ax.text(0.02, 0.95, "Top curated diseases",
            fontsize=10, fontweight="semibold", transform=ax.transAxes, va="top")
    g = dz[dz["gene_symbol"] == gene].nlargest(3, "weighted_sum")
    y_start = 0.80
    max_w = g["weighted_sum"].max() if len(g) else 1.0
    for i, (_, r) in enumerate(g.iterrows()):
        y = y_start - i * 0.25
        name = str(r["disease_name"])
        if len(name) > 40:
            name = name[:38] + "…"
        ax.text(0.02, y + 0.03, name, fontsize=8, color="#333",
                transform=ax.transAxes, va="top")
        bar_len = 0.55 * (r["weighted_sum"] / max_w)
        ax.add_patch(mpatches.Rectangle((0.02, y - 0.08), bar_len, 0.06,
                                         facecolor=GENE_COLOR[gene], alpha=0.85,
                                         transform=ax.transAxes, edgecolor="none"))
        ax.text(0.02 + bar_len + 0.01, y - 0.05,
                f"{r['weighted_sum']:.1f}  (n={int(r['n_evidence'])})",
                fontsize=7.5, color="#111", transform=ax.transAxes, va="top")


def plot(run_dir: Path, out_dir: Path) -> Path:
    apply_style()

    rankings = pd.read_parquet(run_dir / "stage9" / "rankings.parquet")
    verdicts = pd.read_parquet(run_dir / "stage6" / "mechanism_verdicts.parquet")
    gb = pd.read_parquet(run_dir / "stage1" / "1b_genebass" / "genebass_associations.parquet")
    mc = pd.read_parquet(run_dir / "stage3" / "mask_comparison.parquet")
    cv = pd.read_parquet(run_dir / "stage4" / "4a_clinvar" / "clinvar_variants.parquet")
    am = pd.read_parquet(run_dir / "stage2" / "2b_alphamissense" / "alphamissense_scores.parquet")
    hpa = pd.read_parquet(run_dir / "stage8" / "hpa_tissue_ntpm.parquet")
    summary = pd.read_parquet(run_dir / "stage8" / "hpa_expression_summary.parquet")
    dz = pd.read_parquet(run_dir / "stage7" / "curated_per_disease.parquet")

    # 16:9 landscape
    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(3, 6, figure=fig,
                  width_ratios=[1.0, 1.15, 1.1, 1.05, 1.0, 1.15],
                  hspace=0.25, wspace=0.15,
                  left=0.02, right=0.985, top=0.90, bottom=0.04)

    for i, gene in enumerate(GENES):
        rr = rankings[rankings["gene_symbol"] == gene].iloc[0]
        # 0: header
        ax = fig.add_subplot(gs[i, 0])
        _draw_header(ax, gene,
                     verdict=rr["verdict"], rank=int(rr["rank"]),
                     loeuf=float(rr["loeuf"]), composite=float(rr["composite_score"]))
        # 1: Genebass
        ax = fig.add_subplot(gs[i, 1])
        _draw_genebass_card(ax, gb, mc, gene)
        # 2: ClinVar
        ax = fig.add_subplot(gs[i, 2])
        _draw_clinvar_card(ax, cv, gene)
        # 3: AlphaMissense
        ax = fig.add_subplot(gs[i, 3])
        _draw_am_card(ax, am, gene)
        # 4: HPA
        ax = fig.add_subplot(gs[i, 4])
        _draw_hpa_card(ax, hpa, summary, gene)
        # 5: Diseases
        ax = fig.add_subplot(gs[i, 5])
        _draw_disease_card(ax, dz, gene)

    fig.suptitle(
        "Target validation: PCSK9 · JAK2 · LDLR",
        fontsize=16, fontweight="bold", x=0.02, ha="left", y=0.965,
    )
    fig.text(0.02, 0.935,
             "Common-cause target validation on Genebass + ClinVar + AlphaMissense + HPA + Open Targets curated evidence  ·  "
             "PCSK9 and JAK2 → gain-of-function, LDLR → loss-of-function (haploinsufficiency)",
             fontsize=9, color="#555")

    stem = out_dir / "F7_dashboard"
    save(fig, stem)
    plt.close(fig)
    return stem


if __name__ == "__main__":
    run = Path("/mnt/results/rcvtc_runs/sanity_v1")
    outd = run / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    stem = plot(run, outd)
    print(f"wrote {stem}.svg, {stem}.png")
