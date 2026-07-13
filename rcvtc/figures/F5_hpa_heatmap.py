"""F5 — Human Protein Atlas tissue expression heatmap for the 3 targets.

Rows = 3 targets (PCSK9, JAK2, LDLR).
Cols = union of top-8 tissues per gene + core disease-relevant tissues.
Cell = HPA consensus nTPM (log10(1+x) for color scale, raw value labeled in cells).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from figures.palette import GENE_COLOR
from figures.style import apply_style, save

GENES = ["PCSK9", "JAK2", "LDLR"]
MUST_INCLUDE = ["Liver", "Adrenal gland", "Bone marrow", "Blood vessel",
                "Heart muscle", "Kidney", "Lung"]


def _select_tissues(hpa: pd.DataFrame, top_n: int = 8) -> list[str]:
    union: set[str] = set()
    for gene in GENES:
        g = hpa[hpa["gene_symbol"] == gene]
        top = g.nlargest(top_n, "ntpm")["tissue"].tolist()
        union.update(top)
    for t in MUST_INCLUDE:
        if t in hpa["tissue"].values:
            union.add(t)
    # Order tissues by max ntpm across the 3 genes (descending) for a natural gradient
    order = (
        hpa[hpa["tissue"].isin(union)]
        .groupby("tissue")["ntpm"].max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    return order


def plot(run_dir: Path, out_dir: Path) -> Path:
    apply_style()
    hpa = pd.read_parquet(run_dir / "stage8" / "hpa_tissue_ntpm.parquet")
    summary = pd.read_parquet(run_dir / "stage8" / "hpa_expression_summary.parquet")

    tissues = _select_tissues(hpa)

    # Build gene × tissue matrix of ntpm
    mat = np.zeros((len(GENES), len(tissues)))
    for i, gene in enumerate(GENES):
        g = hpa[hpa["gene_symbol"] == gene].set_index("tissue")["ntpm"]
        for j, t in enumerate(tissues):
            mat[i, j] = float(g.get(t, 0.0))

    # log10(1+x) for color scale
    mat_log = np.log10(1.0 + mat)

    fig, ax = plt.subplots(figsize=(11.5, 3.5))

    # Use viridis-r? Prefer a light→dark sequential; use custom (white → gene neutral)
    im = ax.imshow(mat_log, cmap="magma_r", aspect="auto",
                   vmin=0, vmax=max(mat_log.max(), 0.6))

    ax.set_xticks(np.arange(len(tissues)))
    ax.set_xticklabels(tissues, rotation=40, ha="right", fontsize=8.5)
    ax.set_yticks(np.arange(len(GENES)))
    # Colored y-tick labels by gene
    ax.set_yticklabels(GENES, fontsize=10, fontweight="bold")
    for tick, gene in zip(ax.get_yticklabels(), GENES):
        tick.set_color(GENE_COLOR[gene])

    # Annotate cell values (raw nTPM, 1 decimal if <10 else int)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if v < 0.1:
                text = ""
            elif v < 10:
                text = f"{v:.1f}"
            else:
                text = f"{int(round(v))}"
            # Choose text color for contrast
            log_frac = mat_log[i, j] / max(mat_log.max(), 0.6)
            text_color = "white" if log_frac > 0.55 else "#222"
            if text:
                ax.text(j, i, text, ha="center", va="center", fontsize=7, color=text_color)

    ax.set_title("Human Protein Atlas consensus expression (nTPM)", loc="left", pad=6)

    # Add right-side annotation showing tissue specificity summary
    for i, gene in enumerate(GENES):
        s = summary[summary["gene_symbol"] == gene].iloc[0]
        ax.text(
            len(tissues) + 0.4, i,
            f"{s['tissue_specificity']}\n(max: {s['max_tissue']}, {s['max_ntpm']:.1f} nTPM)",
            fontsize=8, va="center", ha="left", color="#333",
        )

    ax.set_xlim(-0.5, len(tissues) + 6.5)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.14, shrink=0.9, aspect=15)
    cbar.set_label("log₁₀(1 + nTPM)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    stem = out_dir / "F5_hpa_heatmap"
    save(fig, stem)
    plt.close(fig)
    return stem


if __name__ == "__main__":
    run = Path("/mnt/results/rcvtc_runs/sanity_v1")
    outd = run / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    stem = plot(run, outd)
    print(f"wrote {stem}.svg, {stem}.png")
