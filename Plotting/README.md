# Plotting

## Active entry point

- `plotting.R`
- `sa_probability_plots.R`

This script regenerates the core four-category figures from the canonical ID files in `Data/`. It expects `Data_wrangling/compile_ssdo_reference_data.R` to have been run first.

The script is now intentionally **self-contained**: the data loading, four-group assignment, summary-table generation, permutation tests, and `ggplot()` code all live directly inside `plotting.R` rather than being split across helper files.

`sa_probability_plots.R` does the same for the three SA probability plots. It now uses the original study-specific inputs in `Data/`: legacy Innocenti p-values, `Assoc_Data.txt` for Ruzicka, and `true_cmh_results.txt` for Neal_Truter_Grieshop, with SNPs mapped onto genes through the local D. melanogaster gene reference.
