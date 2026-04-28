# Data/Summary_tables

This directory keeps only the summary/statistical tables that remain part of the publication workflow.

## Kept table classes

- **Plot Data** tables: precomputed group summaries used to assemble manuscript figures.
- **Permutation Tests** tables: retained statistical test outputs corresponding to the plotted contrasts.
- **Reader-facing analysis tables:** currently `TableS1_SA_Candidate_Enrichment_Tests.csv`, `SA_resolution_enrichment_tests.csv`, `SA_resolution_radar_prevalence.csv`, `SA_resolution_radar_permutation_tests.csv`, and `sa_resolution_feature_notes.csv`.

## Shortcut figure assembly

The kept plot-data tables support direct figure assembly for:

- `Figures/Fig1_SA_Probability_and_Magnitude.png`
- `Figures/Fig2_rmf_tau.png`
- `Figures/Fig3_rlns_tajD_DoS.png`
- `Sup_figures/FigS3_rmf_tau.png`
- `Sup_figures/FigS4_rlns_tajD_DoS.png`

`FigS1`, `FigS2`, and `Fig4` still depend on their active upstream analysis scripts.

Research-only raw-data dumps, overlap bookkeeping tables, and duplicate standalone exports are intentionally no longer kept here.
