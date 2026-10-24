# Plotting and analyses

This directory contains the active main-text plotting scripts together with the manuscript analyses that directly feed those figures and tables.

## Active scripts

1. `boot_permute_metrics.R`  
   Writes the four-category quantitative-metric summaries and permutation tests for *r*<sub>*MF*</sub>, tissue specificity, *pi*<sub>*rN*</sub>, Tajima’s *D*, and *DoS*.

2. `sa_probability_plots.R`  
   Writes the retained Ruzicka SA-probability summary table and within-study permutation tests used in Fig1 and Table1.

3. `sa_magnitude_plots.R`  
   Writes the retained Ruzicka SA-magnitude summary table and within-study permutation tests used in Fig1 and Table1.

4. `figure1_sa_probability_magnitude.R`  
   Builds `Figures/Fig1_SA_Probability_and_Magnitude.png` from the kept study-level summaries.

5. `figure2_figure3_metric_panels.R`  
   Builds:
   - `Figures/Fig2_rmf_tau.png`
   - `Figures/Fig3_rlns_tajD_DoS.png`

6. `sa_candidate_enrichment.R`  
   Runs the FigS1 / TableS1 Fisher enrichment tests for Innocenti and Morrow candidates in their array-assayed universe, Ruzicka et al. candidates in their SNP-mapped universe, and Brigos-Baril et al. (2026) candidates in the protein-coding-gene universe. P-values are BH-adjusted across all nine focal-group-by-candidate-set comparisons.

7. `figure4_sa_resolution.R`  
   Runs the retained study-matched SA-resolution workflow for Fig4 and its downstream summary tables.

8. `figureS8_figure1_with_brigos_5class.R`
   Builds the two-panel FigS8 candidate-class sensitivity figure and its Ruzicka permutation summaries.

The SA probability and magnitude scripts remain separate because they each rebuild a distinct Ruzicka summary dataset and a separate BH-correction family, even though both feed the same manuscript figure.
