"""F2 — Burden β scatter: missense|LC vs pLoF for dual-mask associations.

Every point = one (gene × phenotype) pair with both mask βs.
x-axis: β_missense|LC (per-carrier effect)
y-axis: β_pLoF (per-carrier effect)
size:   -log10(min(p_pLoF, p_missense|LC))
color:  gene
Diagonals: y=x (equal effect), y=1.5x (LoF-haploinsufficiency threshold used in stage 3).

Because both mask βs come from burden averages, missense|LC β is diluted by non-driver
carriers — the y=1.5x rule is a floor, not proof of a LoF mechanism.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from figures.palette import GENE_COLOR, OKABE_ITO
from figures.style import apply_style, save


HERO_LABELS = {
    ("PCSK9", "ldl-direct"): "LDL-C",
    ("LDLR", "date-e78-first-reported-disorders-of-lipoprotein-metabolism-and-other-lipidaemias"): "Lipid disorder",
    ("JAK2", "platelet-count"): "Platelet count",
}


def plot(run_dir: Path, out_dir: Path) -> Path:
    apply_style()
    mc = pd.read_parquet(run_dir / "stage3" / "mask_comparison.parquet")
    d = mc.dropna(subset=["BETA_Burden_missense_LC", "BETA_Burden_pLoF"]).copy()

    # size ∝ -log10(min p)
    minp = np.minimum(d["Pvalue_missense_LC"].values, d["Pvalue_pLoF"].values)
    minp = np.clip(minp, 1e-320, 1.0)
    nlogp = -np.log10(minp)
    # rescale to marker sizes 30-500 for visual range
    size = 30 + 470 * (nlogp - nlogp.min()) / (nlogp.max() - nlogp.min() + 1e-9)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    # scatter per gene for a clean legend
    for gene, sub in d.groupby("gene", sort=False):
        idx = d.index.get_indexer(sub.index)
        ax.scatter(
            sub["BETA_Burden_missense_LC"],
            sub["BETA_Burden_pLoF"],
            s=size[idx],
            color=GENE_COLOR[gene],
            edgecolor="black",
            linewidths=0.7,
            alpha=0.85,
            label=gene,
            zorder=3,
        )

    # Reference diagonals through origin
    x = np.linspace(d["BETA_Burden_missense_LC"].min() * 1.15,
                    d["BETA_Burden_missense_LC"].max() * 1.15, 200)
    ax.plot(x, x, color="black", lw=0.9, ls="--", zorder=1, label="y = x (equal β)")
    ax.plot(x, 1.5 * x, color="black", lw=0.9, ls=":", zorder=1, label="y = 1.5×x (LoF threshold)")

    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)

    # Hero annotations
    for (gene, pheno), label in HERO_LABELS.items():
        row = d[(d["gene"] == gene) & (d["phenotype_slug"] == pheno)]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        # Offset annotation away from the point
        xoff = 0.008 if r["BETA_Burden_missense_LC"] > 0 else -0.008
        yoff = 0.02 if r["BETA_Burden_pLoF"] > 0 else -0.02
        ax.annotate(
            f"{gene} · {label}",
            xy=(r["BETA_Burden_missense_LC"], r["BETA_Burden_pLoF"]),
            xytext=(r["BETA_Burden_missense_LC"] + xoff, r["BETA_Burden_pLoF"] + yoff),
            fontsize=8,
            fontweight="semibold",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.6, shrinkA=0, shrinkB=5),
            zorder=4,
        )

    ax.set_xlabel("β (burden, missense | LC)")
    ax.set_ylabel("β (burden, pLoF)")
    ax.set_title(
        "Genebass burden β per mask\n"
        "Distance below y=x reflects burden-averaging dilution in the missense mask",
        loc="left",
        pad=10,
    )

    # Size legend
    size_marks = [10, 100, 300]
    size_labels = ["10", "100", "300"]
    handles = [
        plt.scatter([], [], s=30 + 470 * (v - nlogp.min()) / (nlogp.max() - nlogp.min() + 1e-9),
                    color="lightgrey", edgecolor="black", linewidths=0.6)
        for v in size_marks
    ]
    size_legend = ax.legend(
        handles, size_labels,
        title="−log10 p (min mask)",
        loc="lower right",
        frameon=True,
        fontsize=8,
        title_fontsize=8,
        borderpad=0.6,
    )
    ax.add_artist(size_legend)

    # Main legend (genes + diagonals)
    ax.legend(loc="upper left", frameon=True, fontsize=9)

    ax.set_aspect("auto")

    stem = out_dir / "F2_burden_scatter"
    save(fig, stem)
    plt.close(fig)
    return stem


if __name__ == "__main__":
    run = Path("/mnt/results/rcvtc_runs/sanity_v1")
    outd = run / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    stem = plot(run, outd)
    print(f"wrote {stem}.svg, {stem}.png")
