# Plotting and analyses

This directory contains the active main-text plotting scripts together with the manuscript analyses that directly feed those figures and tables.

## Active scripts

1. `boot_permute_metrics.R`  
   Writes the four-category quantitative-metric summaries and permutation tests for *r*<sub>*MF*</sub>, tissue specificity, *pi*<sub>*rN*</sub>, Tajima’s *D*, and *DoS*.

2. `sa_magnitude_plots.R`
   Writes the retained Ruzicka SA-effect summary table and within-study permutation tests used in Fig1 and Table1. The metric is the maximum SNP-level Wald chi-squared statistic per protein-coding gene; BH correction covers all six pairwise group contrasts.

3. `figure1_sa_magnitude.R`
   Builds the single-panel `Figures/Fig1_SA_Magnitude.png` from the kept SA-effect summaries.

4. `figure2_figure3_metric_panels.R`
   Builds:
   - `Figures/Fig2_rmf_tau.png`
   - `Figures/Fig3_rlns_tajD_DoS.png`

5. `sa_candidate_enrichment.R`
   Runs the FigS1 / TableS1 Fisher enrichment tests for Innocenti and Morrow candidates in their array-assayed universe, Ruzicka et al. candidates in their SNP-mapped universe, and Brigos-Baril et al. (2026) candidates in the protein-coding-gene universe. TableS1 retains raw two-sided P-values and adds BH-adjusted P-values across all nine tests. FigS1 asterisks indicate BH-adjusted P < 0.05 across that same nine-test family.

   Each bar shows the percentage of study candidates within the focal group or its complement in that column's universe, not the percentage of candidates belonging to the focal group. For example, the Innocenti and Morrow SSDO bar is 41/560 (7.3%); 560 is the number of array-assayed SSDOs, while the matched candidate set contains 915 genes. Disease orthologues exclude SSDOs, and orthologues exclude disease orthologues; "other genes" is the complement of the focal group, not just the Background category.

6. `figure4_sa_resolution.R`
   Runs the retained study-matched SA-resolution workflow for Fig4 and its downstream summary tables.

7. `figureS4_figure1_with_brigos_5class.R`
   Builds the single-panel `Sup_figures/FigS4_SA_Magnitude_Brigos5.png` candidate-class sensitivity figure and its Ruzicka SA-effect permutation summaries.
