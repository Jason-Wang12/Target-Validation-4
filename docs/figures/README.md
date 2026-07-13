# rcvtc figures — sanity_v1

Publication-quality figures generated from the `sanity_v1` run of the
`rcvtc` common-cause target-validation pipeline. Every value shown in
these figures is read directly from the parquet outputs under
`rcvtc_runs/sanity_v1/`. No numbers are hard-coded.

## Figures

### F2 — Burden β scatter (`F2_burden_scatter.{svg,png}`)
Genebass burden β for the pLoF mask vs. the missense|LC mask, one point
per (gene × phenotype) association where both masks have signal
(n=25 dual-mask associations). Point size scales with -log10 of the more
significant p-value. Reference diagonals: y=x (equal β) and y=1.5x
(loss-of-function-haploinsufficiency threshold used by the stage-3
classifier). Hero labels on PCSK9 LDL-C, LDLR lipid disorder,
and JAK2 platelet count.

**Interpretation:** Every dual-mask signal falls above y=1.5x, meaning
pLoF β magnitude always exceeds missense|LC β magnitude in Genebass
burden output — a burden-averaging artifact where a single missense
driver (e.g. JAK2 V617F, PCSK9 D374Y) is diluted by non-driver missense
carriers. This is why stage 3 alone cannot distinguish GoF from LoF;
stage 6 combines ClinVar spectrum + LOEUF + curated GoF flags to make
the final call.

*Data sources:* `stage3/mask_comparison.parquet`.

### F4 — AlphaMissense landscape (`F4_am_landscape.{svg,png}`)
Three vertically stacked panels (one per gene). Each panel shows:
- **Top strip:** InterPro domain/family regions (colored rectangles).
- **Main plot:** AlphaMissense per-residue max pathogenicity as vertical
  lines colored by AM class (grey = likely_benign, yellow = ambiguous,
  vermilion = likely_pathogenic). Reference lines at 0.34 and 0.564
  (AM class boundaries).
- **Bottom strip:** ClinVar 2⁺★ pathogenic variant positions (blue =
  missense, black = LoF, green = other/unknown), with count summary.

**Interpretation:** Compare the per-gene ClinVar footprint to the AM
predictions. PCSK9 pathogenic missense variants (E32K, D374Y/H, R215H)
cluster around the catalytic domain (156-421); LDLR pathogenic LoF
alleles (n=86) are scattered but concentrated in the LDLa repeats and
EGF-like regions; JAK2 shows a single 2⁺★ pathogenic variant (V617F,
position ~617 in the JH2 pseudokinase domain).

*Data sources:* `stage2/2b_alphamissense/alphamissense_scores.parquet`,
`stage4/4a_clinvar/clinvar_variants.parquet`,
`stage2/2c_uniprot/interpro_domains_compact.parquet`.

### F5 — HPA tissue heatmap (`F5_hpa_heatmap.{svg,png}`)
3 × 20 heatmap of HPA consensus nTPM (log10(1+x) color scale, raw values
labeled in cells). Tissue set is the union of the top-8 tissues per gene
plus 7 disease-relevant tissues (Liver, Adrenal gland, Bone marrow,
Blood vessel, Heart muscle, Kidney, Lung). Right-hand summary shows the
HPA tissue-specificity call and max tissue per gene.

**Interpretation:**
- PCSK9 (Tissue enriched, liver): Liver 47.5 nTPM (matches secreted
  hepatic biology; primary target of alirocumab/evolocumab).
- JAK2 (Low tissue specificity): highest in blood vessel / heart
  muscle / bone marrow, reflecting broad hematopoietic and vascular
  cytokine signaling.
- LDLR (Tissue enhanced, adrenal gland): Adrenal 112.1, Liver 56.4;
  LDL-C clearance biology is dominated by hepatic expression despite the
  higher absolute value in adrenal.

*Data sources:* `stage8/hpa_tissue_ntpm.parquet`,
`stage8/hpa_expression_summary.parquet`.

### F6 — Top curated diseases (`F6_top_diseases.{svg,png}`)
Three horizontal bar panels showing the top-6 Open Targets curated
disease associations per gene, ranked by `weighted_sum` (sum of
per-evidence scores across curated data sources). Bar labels show
`weighted_sum` and `n_evidence`.

**Interpretation:**
- PCSK9 and LDLR share the familial hypercholesterolemia / lipid
  metabolism disease axis (LDLR carries stronger evidence: weighted_sum
  ~190 vs. ~107).
- JAK2 diseases (myeloproliferative disorder, polycythemia, familial
  thrombocytosis) reflect the V617F GoF gain-of-signaling mechanism.

*Data sources:* `stage7/curated_per_disease.parquet`.

### F7 — Target validation dashboard (`F7_dashboard.{svg,png}`)
Landscape 16:9 composite. 3 gene rows × 6 evidence columns (Header /
Genebass burden / ClinVar / AlphaMissense / HPA / Diseases). Header
shows final verdict pill, LOEUF, and composite score from stage 9.
All values are pulled from the source parquets.

**Interpretation:** Reading a row left-to-right gives the full evidence
stack behind each verdict:
- LDLR (rank 1, LoF · haploinsufficiency): Genebass shows large pLoF
  β (+0.14) alongside 86 LoF alleles in ClinVar and 359 pathogenic
  missense variants — consistent with classical loss-of-function.
- PCSK9 (rank 2, GoF · activating): Genebass burden shows pLoF β
  larger than missense β, but all 5 ClinVar pathogenic variants are
  missense (0 LoF) and LOEUF=1.14 indicates LoF tolerance — pointing to
  a GoF mechanism driven by specific residues (D374Y, R215H, S127R).
- JAK2 (rank 3, GoF · activating): 1 ClinVar pathogenic (V617F),
  0 LoF, curated GoF term present, LOEUF=0.74 — the burden signal is
  averaged across non-driver carriers and diluted, but the mechanism is
  the canonical V617F pseudokinase activating mutation.

*Data sources:* `stage9/rankings.parquet`,
`stage6/mechanism_verdicts.parquet`, `stage1/1b_genebass/`,
`stage3/mask_comparison.parquet`, `stage4/4a_clinvar/`,
`stage2/2b_alphamissense/`, `stage8/hpa_*.parquet`,
`stage7/curated_per_disease.parquet`.

## Notes on data and figures

- All figures use the Okabe-Ito colorblind-friendly palette
  (`PCSK9=#E69F00`, `LDLR=#0072B2`, `JAK2=#CC79A7`).
- Files are exported both as SVG (editable text, `svg.fonttype='none'`)
  and PNG (300 DPI). Fonts: Liberation Sans / DejaVu Sans.
- Figure code lives in `/mnt/results/rcvtc/figures/`
  (`palette.py`, `style.py`, `interpro.py`, and one file per figure).

## Reproducibility

```bash
cd /mnt/results/rcvtc
python figures/F2_burden_scatter.py
python figures/F4_am_landscape.py
python figures/F5_hpa_heatmap.py
python figures/F6_top_diseases.py
python figures/F7_dashboard.py
```
