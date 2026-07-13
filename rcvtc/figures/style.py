"""Matplotlib style: Liberation Sans + editable SVG text + tidy defaults."""
import matplotlib as mpl

def apply_style():
    mpl.rcParams.update({
        "font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.linewidth": 0.9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "svg.fonttype": "none",       # editable SVG text
        "pdf.fonttype": 42,           # embed TrueType
        "ps.fonttype": 42,
    })


def save(fig, path_stem, dpi=300):
    """Save figure as SVG + PNG at path_stem (no extension)."""
    from pathlib import Path
    p = Path(path_stem)
    p.parent.mkdir(parents=True, exist_ok=True)
    svg_path = p.with_suffix(".svg")
    png_path = p.with_suffix(".png")
    fig.savefig(svg_path, bbox_inches="tight", dpi=dpi)
    fig.savefig(png_path, bbox_inches="tight", dpi=dpi)
    return {"svg": str(svg_path), "png": str(png_path)}
