# Data/Summary_tables

This directory keeps only the summary/statistical tables that remain part of the publication workflow.

## Kept table classes

- **Plot Data** tables: precomputed group summaries used to assemble manuscript figures.
- **Permutation Tests** tables: retained statistical test outputs corresponding to the plotted contrasts.
- **Reader-facing analysis tables:** currently `TableS1_SA_Candidate_Enrichment_Tests.csv`, `TableS1_SA_Candidate_Enrichment_Overlaps.csv`, `SA_resolution_enrichment_tests.csv`, `SA_resolution_multivariate_tests.csv`, `SA_resolution_radar_prevalence.csv`, `SA_resolution_radar_permutation_tests.csv`, `SA_resolution_tse_thresholds.csv`, and `sa_resolution_feature_notes.csv`. The candidate-enrichment tests include all three candidate sets, raw two-sided Fisher `PValue`, and `PAdjustedBH` across all nine focal-group by candidate-study tests; `AdjustmentFamily` records that family. FigS1 asterisks use `PAdjustedBH < 0.05`. The SA-resolution tables use study-matched denominators and record `NConsidered` where applicable.

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

The `FigS7_rmf_*` and `FigS7_sa_effect_*` exports contain gene-inclusion audits, stage-specific gene values, bootstrap summaries, and permutation tests for FigS7. Penultimate is MB or TS; Ultimate is HMB or HTS. Categories are mutually exclusive, but stages may overlap; shared genes remain in both stages. Each metric's `permutation_tests.csv` has four within-category two-sided comparisons with 10,000 permutations, holding shared genes fixed and shuffling exclusive memberships. BH correction is separate within rMF and within SA effects, across four tests each. `NGroup1`, `NGroup2`, `NOverlap`, exclusive counts, and `NUniqueGenes` record the tested universes. `FigS7_genome_wide_reference_means.csv` records the grey reference lines: arithmetic means across all protein-coding genes with finite metric values, not just genes in either stage or the intersection of metric assays. The former FigS7 radar-overlay test is archived with its script and figure in `side/archived_figS7_radar_context/`.

`FigS8_category_coverage.csv` records canonical eligible, excluded-before-alignment, and measured genes by category, plus one-to-one measured counts and, when provided by the refreshed mapping output, one-to-one eligibility and missingness. `FigS8_group_bootstrap_summaries.csv` contains means and percentile bootstrap intervals for whole-alignment identity and fly/human coverage in the main and one-to-one sets. `FigS8_pairwise_permutation_tests.csv` contains the six within-metric BH families (three category contrasts per metric and analysis set), including sample sizes and shared human partners. These are exploratory conservation analyses; the current exports reflect the initial mapping filter until the pending reference-aware recovery is completed.

Research-only raw-data dumps, overlap bookkeeping tables, and duplicate standalone exports are intentionally no longer kept here.
