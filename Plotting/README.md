# Main-figure plotting

This directory contains the active scripts for the retained **main-figure** plotting workflow.

## Workflow split

- `plotting.R`, `sa_probability_plots.R`, and `sa_magnitude_plots.R` generate the kept **plot-data** and **permutation-test** tables in `Data/Summary_tables/`.
- `figure1_sa_probability_magnitude.R` and `figure2_figure3_metric_panels.R` assemble the publication figures in `Figures/`.

## Active scripts

1. `plotting.R`  
   Writes the four-category metric summaries and permutation tests for *r*<sub>*MF*</sub>, tissue specificity, *pi*<sub>*rN*</sub>, Tajima’s *D*, and *DoS*.

2. `sa_probability_plots.R`  
   Writes the kept Ruzicka and Neal/Truter/Grieshop SA-probability summary tables.

3. `sa_magnitude_plots.R`  
   Writes the kept Ruzicka and Neal/Truter/Grieshop SA-magnitude summary tables.

4. `figure1_sa_probability_magnitude.R`  
   Builds `Figures/Fig1_SA_Probability_and_Magnitude.png`. If the kept summary tables already exist, this script can reuse them directly instead of recomputing the study summaries.

5. `figure2_figure3_metric_panels.R`  
   Builds:
   - `Figures/Fig2_rmf_tau.png`
   - `Figures/Fig3_rlns_tajD_DoS.png`

Standalone component PNGs are no longer part of the tracked publication surface.
