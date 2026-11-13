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

`FigS1`, `FigS2`, `FigS3`, `Fig4`, `FigS4`, and `FigS9` still depend on their active upstream analysis scripts.

The `Ruzicka_et_al_SA_magnitude_*` tables use `max(abs(Wald_Stat))` over finite SNP statistics overlapping each protein-coding gene. Absolute values are taken before the gene maximum, not afterwards; no squaring is applied. This maximum absolute Wald z-statistic is the SA-effect metric throughout the active analyses. The four-class tests support Table1; the `_brigos5` tables support FigS4. `TableS2_candidate_context_permutation_tests.csv` contains the three SSDO-versus-candidate SA-effect comparisons for FigS3 (Innocenti and Morrow, Grieshop, and Ruzicka), using overlap-aware permutation tests and BH correction across those three comparisons.

The `FigS3_*.csv` exports retain plotted gene-level values, four-group bootstrap summaries, reference-line values, and a Grieshop candidate inclusion audit. The audit records all 695 significant DE genes, protein-coding membership, mapped Ruzicka values, and inclusion status; 530 are plotted. Gene-level rows are unique within each candidate group but can occur in multiple groups.

The `FigS9_rmf_*` and `FigS9_sa_effect_*` exports contain gene-inclusion audits, stage-specific gene values, bootstrap summaries, and permutation tests for FigS9. Penultimate is MB or TS; Ultimate is HMB or HTS. Categories are mutually exclusive, but stages may overlap; shared genes remain in both stages. Each metric's `permutation_tests.csv` has four within-category two-sided comparisons with 10,000 permutations, holding shared genes fixed and shuffling exclusive memberships. BH correction is separate within rMF and within SA effects, across four tests each. `NGroup1`, `NGroup2`, `NOverlap`, exclusive counts, and `NUniqueGenes` record the tested universes. `FigS9_genome_wide_reference_means.csv` records the grey reference lines: arithmetic means across all protein-coding genes with finite metric values, not just genes in either stage or the intersection of metric assays. The former FigS7 radar-overlay test is archived with its script and figure in `side/archived_figS7_radar_context/`.

`FigS7_category_coverage.csv` records canonical eligible, excluded-before-alignment, and measured genes by category. `FigS7_group_bootstrap_summaries.csv` contains means and percentile bootstrap intervals for whole-alignment identity and fly/human coverage in all measured orthologues. `FigS7_pairwise_permutation_tests.csv` supplies TableS6 and contains three within-metric BH families (three category contrasts per metric), including sample sizes and shared human partners. No one-to-one filtered analyses are performed. These remain exploratory conservation analyses; 29 of 9,092 eligible genes lack a validated protein pair and are not imputed.

`FigS8_*.csv` records conservation-adjusted SA regression inputs and results.
`gene_inclusion_audit` lists all protein-coding genes and explicit exclusions;
`category_coverage` separates SA and conservation missingness;
`gene_level_values` retains fitted values and residuals;
`summaries` contains raw/adjusted means with human-partner-cluster bootstrap
intervals; `coefficients` and `diagnostics` record the model, covariate
centering/scaling, condition numbers, residual variability, and R-squared.
`tests` compares global three-category and disease-only-versus-other effects
before and after adjustment, using restricted wild-cluster null resampling
(10,000 replicates). BH covers the four tests in the full measured SA/conservation universe.
`input_manifest` and `FigS8_R_session.txt` record data/script hashes and the R
environment. Background is not assigned a conservation value. All raw versus
adjusted comparisons use identical SA/conservation-matched genes.

`SA_effect_max_abs_revision_audit.csv` is the fixed 2026-11-12 revision audit:
all 27 SA contrasts before (`max(Wald_Stat)`) and after
(`max(abs(Wald_Stat))`) the correction, with observed statistics, raw and
BH-adjusted P-values, sample sizes/overlap, adjustment families, and S9
contrast intervals. Gene universes, resampling methods, and seeds were held
fixed. Only the Fig1/Table1 orthologues-versus-background comparison changed
its BH < 0.05 decision. This historical comparison is not regenerated by the
figure scripts.

Research-only raw-data dumps, overlap bookkeeping tables, and duplicate standalone exports are intentionally no longer kept here.
