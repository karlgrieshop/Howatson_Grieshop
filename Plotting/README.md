# Plotting

## Active entry point

- `plotting.R`
- `sa_probability_plots.R`
- `sa_magnitude_plots.R`

These scripts regenerate the active figures from the canonical ID files in `Data/outputs/` and the raw metric inputs in `Data/inputs/`. They expect `Data_wrangling/compile_ssdo_reference_data.R` to have been run first.

The script is now intentionally **self-contained**: the data loading, four-group assignment, summary-table generation, permutation tests, and `ggplot()` code all live directly inside `plotting.R` rather than being split across helper files.

`sa_probability_plots.R` now produces the default protein-coding-universe SA probability plots in `Figures/`, plus full-universe sensitivity versions in `Sup_figures/`, for the Ruzicka and Neal_Truter_Grieshop datasets.

`sa_magnitude_plots.R` produces the analogous SA magnitude plots and permutation summaries, using `Wald_Stat` for Ruzicka SNPs and finite `abs(LogOdds)` for Neal_Truter_Grieshop SNPs before collapsing to genes.
