"""Okabe-Ito colorblind-friendly palette for rcvtc figures."""

# Okabe-Ito 8-color palette
OKABE_ITO = {
    "black":   "#000000",
    "orange":  "#E69F00",
    "sky":     "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermilion": "#D55E00",
    "purple":  "#CC79A7",
    "grey":    "#999999",
}

# Gene → color (consistent across all figures)
GENE_COLOR = {
    "PCSK9": OKABE_ITO["orange"],
    "LDLR":  OKABE_ITO["blue"],
    "JAK2":  OKABE_ITO["purple"],
}

# AlphaMissense class colors
AM_CLASS_COLOR = {
    "likely_benign":     OKABE_ITO["grey"],
    "ambiguous":         OKABE_ITO["yellow"],
    "likely_pathogenic": OKABE_ITO["vermilion"],
}

# Verdict color
VERDICT_COLOR = {
    "lof_haploinsufficiency": OKABE_ITO["blue"],
    "gof_activating":         OKABE_ITO["vermilion"],
    "gof_candidate":          OKABE_ITO["orange"],
    "missense_dominant":      OKABE_ITO["purple"],
    "lof_probable":           OKABE_ITO["sky"],
    "undetermined":           OKABE_ITO["grey"],
    "single_mask_only":       OKABE_ITO["grey"],
}
