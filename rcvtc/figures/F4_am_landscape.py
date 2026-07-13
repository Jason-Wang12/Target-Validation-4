"""F4 — AlphaMissense pathogenicity landscape with ClinVar + InterPro tracks.

Three vertically stacked panels (one per gene: PCSK9, JAK2, LDLR).

Per panel:
  Top strip:    InterPro domain rectangles (labeled inside where they fit)
  Main:         AM per-residue max pathogenicity (grey/yellow/red by am_class)
  Bottom strip: ClinVar 2+ star pathogenic ticks (missense=blue, LoF=red, unknown=grey)

x-axis: protein position (residue)
y-axis (main): AM pathogenicity ∈ [0, 1]
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
from figures.palette import AM_CLASS_COLOR, GENE_COLOR
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


def _extract_position(pc) -> int | None:
    if pc is None or (isinstance(pc, float) and np.isnan(pc)):
        return None
    tok = str(pc).split(",")[0].strip()
    m = re.search(r"[A-Z][a-z]{0,2}(\d+)", tok)
    if m:
        return int(m.group(1))
    return None


CLASSIFY_COLOR = {
    "missense": "#0072B2",  # Okabe-Ito blue
    "lof":      "#000000",  # black — distinct from AM likely_pathogenic (vermilion)
    "unknown":  "#009E73",  # Okabe-Ito bluish green — distinct from AM likely_benign (grey)
}


def _draw_panel(fig, gs_outer, gene: str, am_g: pd.DataFrame, cv_g: pd.DataFrame, ip_g: pd.DataFrame):
    """Build a 3-row nested gridspec for one gene."""
    gs = gs_outer.subgridspec(3, 1, height_ratios=[1.1, 4.0, 0.9], hspace=0.05)
    ax_dom = fig.add_subplot(gs[0, 0])
    ax_am = fig.add_subplot(gs[1, 0], sharex=ax_dom)
    ax_cv = fig.add_subplot(gs[2, 0], sharex=ax_dom)

    # -------- InterPro domain track --------
    ax_dom.set_ylim(0, 2.4)
    ax_dom.set_yticks([])
    ax_dom.set_facecolor("white")
    for spine in ("top", "left", "right"):
        ax_dom.spines[spine].set_visible(False)
    ax_dom.spines["bottom"].set_visible(False)
    ax_dom.tick_params(bottom=False, labelbottom=False)

    dom_color = GENE_COLOR[gene]
    for _, row in ip_g.iterrows():
        rect = mpatches.Rectangle(
            (row["start"], 0.15), row["end"] - row["start"], 0.55,
            facecolor=dom_color, alpha=0.55, edgecolor="black", linewidth=0.5,
        )
        ax_dom.add_patch(rect)
        # Label wide boxes inside; small boxes get a label just above at 45°
        width = row["end"] - row["start"]
        if width > 40:
            ax_dom.text((row["start"] + row["end"]) / 2, 0.42, row["short_label"],
                        ha="center", va="center", fontsize=6.5, fontweight="semibold")
        else:
            ax_dom.text((row["start"] + row["end"]) / 2, 0.85, row["short_label"],
                        ha="left", va="bottom", fontsize=5.5, color="#333",
                        rotation=45, rotation_mode="anchor")

    # -------- AM main --------
    # Aggregate to per-residue max pathogenicity + assign class by max-score entry
    if len(am_g):
        am_pos = am_g.copy()
        am_pos = am_pos.dropna(subset=["protein_start", "am_pathogenicity"])
        # per-residue max
        agg = am_pos.sort_values("am_pathogenicity", ascending=False).drop_duplicates("protein_start")
        # Plot as vertical lines by class for visual density
        for cls, color in AM_CLASS_COLOR.items():
            sub = agg[agg["am_class"] == cls]
            if len(sub) == 0:
                continue
            ax_am.vlines(sub["protein_start"], 0, sub["am_pathogenicity"],
                         color=color, alpha=0.7, linewidth=0.7)
        max_pos = int(am_pos["protein_start"].max())
    else:
        max_pos = int(ip_g["end"].max()) if len(ip_g) else 500

    ax_am.set_ylim(0, 1.05)
    ax_am.set_ylabel("AM pathogenicity", fontsize=9)
    ax_am.axhline(0.564, color="grey", lw=0.5, ls=":", zorder=0)  # AM likely-pathogenic threshold
    ax_am.axhline(0.34, color="grey", lw=0.5, ls=":", zorder=0)   # ambiguous cutoff
    ax_am.tick_params(labelbottom=False)
    ax_am.set_facecolor("white")
    for spine in ("top", "right"):
        ax_am.spines[spine].set_visible(False)

    # Gene label + count of likely-pathogenic
    n_path = int((am_g["am_class"] == "likely_pathogenic").sum()) if len(am_g) else 0
    ax_am.text(0.005, 0.95, f"{gene}",
               transform=ax_am.transAxes, fontsize=11, fontweight="bold", va="top", color=GENE_COLOR[gene])
    ax_am.text(0.005, 0.88, f"AM likely-pathogenic residues: {n_path}",
               transform=ax_am.transAxes, fontsize=8, va="top", color="#333")

    # -------- ClinVar bottom --------
    ax_cv.set_ylim(0, 1)
    ax_cv.set_yticks([])
    for spine in ("top", "left", "right"):
        ax_cv.spines[spine].set_visible(False)
    ax_cv.set_facecolor("white")

    path_g = cv_g[cv_g["clinsig_bucket"] == "pathogenic_lp"].copy()
    n_path_missense = 0
    n_path_lof = 0
    n_path_unknown = 0
    for _, r in path_g.iterrows():
        pos = _extract_position(r["protein_change"])
        if pos is None:
            continue
        cls = _classify_protein_change(r["protein_change"])
        if cls == "missense":
            n_path_missense += 1
        elif cls == "lof":
            n_path_lof += 1
        else:
            n_path_unknown += 1
        ax_cv.vlines(pos, 0.15, 0.85, color=CLASSIFY_COLOR[cls], linewidth=0.7, alpha=0.85)

    ax_cv.set_xlim(0, max_pos + 10)
    ax_cv.set_xlabel("Protein position (residue)", fontsize=9)
    # Right-anchor the summary text so it doesn't collide with dense N-terminal ticks (LDLR)
    ax_cv.text(0.995, 1.35,
               f"ClinVar 2⁺★ pathogenic: n={len(path_g)}  (missense={n_path_missense}, LoF={n_path_lof}, other={n_path_unknown})",
               transform=ax_cv.transAxes, fontsize=8, ha="right", va="top", color="#333")

    return ax_dom, ax_am, ax_cv


def plot(run_dir: Path, out_dir: Path) -> Path:
    apply_style()

    am = pd.read_parquet(run_dir / "stage2" / "2b_alphamissense" / "alphamissense_scores.parquet")
    cv = pd.read_parquet(run_dir / "stage4" / "4a_clinvar" / "clinvar_variants.parquet")
    ip = pd.read_parquet(run_dir / "stage2" / "2c_uniprot" / "interpro_domains_compact.parquet")

    fig = plt.figure(figsize=(11.5, 9.5))
    gs = GridSpec(3, 1, figure=fig, hspace=0.45)

    for i, gene in enumerate(GENES):
        _draw_panel(
            fig, gs[i, 0], gene,
            am_g=am[am["gene_symbol"] == gene],
            cv_g=cv[cv["gene_symbol"] == gene],
            ip_g=ip[ip["gene_symbol"] == gene],
        )

    # Legends
    legend_ax = fig.add_axes([0.15, 0.965, 0.75, 0.02])
    legend_ax.axis("off")
    handles_am = [mpatches.Patch(color=c, label=lbl.replace("_", " "))
                  for lbl, c in AM_CLASS_COLOR.items()]
    handles_cv = [mpatches.Patch(color=CLASSIFY_COLOR["missense"], label="ClinVar missense"),
                  mpatches.Patch(color=CLASSIFY_COLOR["lof"], label="ClinVar LoF"),
                  mpatches.Patch(color=CLASSIFY_COLOR["unknown"], label="ClinVar other")]
    legend_ax.legend(
        handles=handles_am + handles_cv,
        loc="center",
        ncol=6,
        frameon=False,
        fontsize=8,
    )

    fig.suptitle("AlphaMissense per-residue pathogenicity, ClinVar 2⁺★ pathogenic variants, and InterPro domains",
                 fontsize=11, x=0.02, ha="left", y=0.995, weight="semibold")

    stem = out_dir / "F4_am_landscape"
    save(fig, stem)
    plt.close(fig)
    return stem


if __name__ == "__main__":
    run = Path("/mnt/results/rcvtc_runs/sanity_v1")
    outd = run / "figures"
    outd.mkdir(parents=True, exist_ok=True)
    stem = plot(run, outd)
    print(f"wrote {stem}.svg, {stem}.png")
