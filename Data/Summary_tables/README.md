# Data/Summary_tables

This directory keeps only the summary/statistical tables that remain part of the publication workflow.

## Kept table classes

- **Plot Data** tables: precomputed group summaries used to assemble manuscript figures.
- **Permutation Tests** tables: retained statistical test outputs corresponding to the plotted contrasts.
- **Reader-facing analysis tables:** currently `TableS1_SA_Candidate_Enrichment_Tests.csv`, `TableS1_SA_Candidate_Enrichment_Overlaps.csv`, `SA_resolution_enrichment_tests.csv`, `SA_resolution_multivariate_tests.csv`, `SA_resolution_radar_prevalence.csv`, `SA_resolution_radar_permutation_tests.csv`, `SA_resolution_tse_thresholds.csv`, and `sa_resolution_feature_notes.csv`. The candidate-enrichment tests include all three candidate sets and nominal two-sided Fisher P-values without multiple-testing adjustment, matching the descriptive annotations in FigS1/TableS1; the SA-resolution tables use study-matched denominators and record `NConsidered` where applicable.

## Shortcut figure assembly

The kept plot-data tables support direct figure assembly for:

- `Figures/Fig1_SA_Magnitude.png`
- `Figures/Fig2_rmf_tau.png`
- `Figures/Fig3_rlns_tajD_DoS.png`
- `Sup_figures/FigS5_rmf_tau.png`
- `Sup_figures/FigS6_rlns_tajD_DoS.png`

`FigS1`, `FigS2`, `FigS3`, `Fig4`, `FigS4`, and `FigS7` still depend on their active upstream analysis scripts.

The `Ruzicka_et_al_SA_magnitude_*` tables retain the per-gene maximum Wald chi-squared metric used for SA-effect comparisons. The four-class tests support Table1; the `_brigos5` tables support FigS4. `TableS2_candidate_context_permutation_tests.csv` contains the three SSDO-versus-candidate SA-effect comparisons for FigS3 (Innocenti and Morrow, Grieshop, and Ruzicka), using overlap-aware permutation tests and BH correction across those three comparisons.

The `FigS3_*.csv` exports retain plotted gene-level values, four-group bootstrap summaries, reference-line values, and a Grieshop candidate inclusion audit. The audit records all 695 significant DE genes, protein-coding membership, mapped Ruzicka values, and inclusion status; 530 are plotted. Gene-level rows are unique within each candidate group but can occur in multiple groups.

Research-only raw-data dumps, overlap bookkeeping tables, and duplicate standalone exports are intentionally no longer kept here.
