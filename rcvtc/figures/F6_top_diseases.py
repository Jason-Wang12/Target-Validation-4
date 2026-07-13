"""F6 — Top curated disease associations per target (Open Targets curated).

Three horizontal bar panels (one per gene). Bar length = weighted_sum
(sum of per-evidence scores across data sources). Labels show n_evidence.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from figures.palette import GENE_COLOR
from figures.style import apply_style, save

GENES = ["PCSK9", "JAK2", "LDLR"]
TOP_N = 6


def _wrap_disease(name: str, width: int = 40) -> str:
    """Wrap long disease names at word boundaries."""
    words = str(name).split()
    lines = []
    cur = ""
    for w in words:
        if len(cur) + len(w) + 1 <= width:
            cur = (cur + " " + w).strip()
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return "\n".join(lines)


def plot(run_dir: Path, out_dir: Path) -> Path:
    apply_style()
    d = pd.read_parquet(run_dir / "stage7" / "curated_per_disease.parquet")

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(9.5, 8.0),
                             gridspec_kw={"hspace": 0.55})

    for ax, gene in zip(axes, GENES):
        top = d[d["gene_symbol"] == gene].nlargest(TOP_N, "weighted_sum")
        top = top.iloc[::-1]  # so largest at top after barh flips

        labels = [_wrap_disease(x, width=40) for x in top["disease_name"]]
        y = range(len(top))
        color = GENE_COLOR[gene]

        bars = ax.barh(y, top["weighted_sum"], color=color, alpha=0.85, edgecolor="black", linewidth=0.4)

        for i, (bar, (_, r)) in enumerate(zip(bars, top.iterrows())):
            w = bar.get_width()
            ax.text(w + max(top["weighted_sum"]) * 0.015, i,
                    f"{w:.1f}  (n={int(r['n_evidence'])})",
                    va="center", fontsize=8, color="#222")

        ax.set_yticks(list(y))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0, max(top["weighted_sum"]) * 1.30)
        ax.set_title(gene, loc="left", fontsize=11, fontweight="bold",
                     color=GENE_COLOR[gene], pad=4)
        ax.set_xlabel("weighted_sum (Σ per-evidence score across curated sources)", fontsize=8)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle("Top curated disease associations per target (Open Targets curated evidence)",
                 fontsize=11, weight="semibold", x=0.02, ha="left", y=0.995)

    stem = out_dir / "F6_top_diseases"
    save(fig, stem)
    plt.close(fig)
    return stem


if __name__ == "__main__":
    run = Path("/mnt/results/rcvtc_runs/sanity_v1")
    outd = run / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    stem = plot(run, outd)
    print(f"wrote {stem}.svg, {stem}.png")
