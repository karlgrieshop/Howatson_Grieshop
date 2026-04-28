# Supplemental plotting

- `sbge_metric_plots.R` regenerates the supplemental four-category metric figures using a single SBGE x-axis with Mishra whole-body sex-bias classes plus `All` at the far left.
- `figureS3_figureS4_metric_panels.R` packages those SBGE summary tables into the stacked supplemental composites `FigS3_rmf_tau.png` and `FigS4_rlns_tajD_DoS.png`.
- It uses the current canonical SSDO / disease-orthologue / orthologue ID sets and the current bootstrap + permutation workflow, then subsets to genes present in `Data/outputs/mishra_whole_body_sbge_bins.csv`.
- The six x-axis bins are:
  - All
  - Highly female-biased
  - Female-biased
  - Unbiased
  - Male-biased
  - Highly male-biased

`Female-biased` and `Male-biased` in these supplemental panels correspond to Mishra's `Moderately female-biased` and `Moderately male-biased` categories from the rebuilt MaleLimitedEvo-style binning.
