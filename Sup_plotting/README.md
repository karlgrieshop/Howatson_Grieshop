# Supplemental plotting

- `sbge_metric_plots.R` regenerates the supplemental four-category metric figures using a single SBGE x-axis with Mishra whole-body sex-bias classes plus `All` at the far left.
- It uses the current canonical SSDO / disease-orthologue / orthologue ID sets and the current bootstrap + permutation workflow, then subsets to genes present in `Data/Sex.Biased.Genes.Whole.csv`.
- The six x-axis bins are:
  - All
  - Highly female-biased
  - Female-biased
  - Unbiased
  - Male-biased
  - Highly male-biased

`Female-biased` and `Male-biased` in these supplemental panels correspond to Mishra's `Moderately female-biased` and `Moderately male-biased` categories.
