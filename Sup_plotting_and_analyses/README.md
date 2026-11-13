# Supplementary plotting and analyses

This directory contains the active scripts for the retained **supplementary figure** plotting workflow.

## Active scripts

1. `sbge_metric_plots.R`  
   Writes the SBGE-stratified metric plot-data and permutation-test tables kept in `Data/Summary_tables/SBGE_supplement/`.

2. `figureS2_sa_mag.R`
   Builds the single-panel `Sup_figures/FigS2_sa_mag.png` SA-effect comparison stratified by whole-body sex-biased gene expression.

3. `figureS3_context.R`
   Builds the single-panel `Sup_figures/FigS3_sa_mag_context.png`, comparing SSDOs with Innocenti and Morrow (2010), Grieshop et al. (2025), and Ruzicka et al. (2019) candidates using per-gene maximum absolute Ruzicka Wald z-statistics. Grieshop candidates are the 695 `Sig == TRUE` genes in `Data/inputs/DE.summary.tsv`; 530 have protein-coding, SNP-mapped estimates. Candidate groups may overlap. Horizontal blue/green/grey reference lines show Fig1-category bootstrap means for disease orthologues, orthologues, and background, using the same corrected metric. The script exports `FigS3_gene_level_values.csv`, `FigS3_bootstrap_summaries.csv`, `FigS3_reference_lines.csv`, and `FigS3_Grieshop_candidate_audit.csv` to `Data/Summary_tables/`.

4. `figureS5_figureS6_metric_panels.R`
   Builds:
   - `Sup_figures/FigS5_rmf_tau.png`
   - `Sup_figures/FigS6_rlns_tajD_DoS.png`

5. `figureS7_resolution_stages.R`
   Builds `Sup_figures/FigS7_resolution_stages.png`: rMF (A) and SA effects (B) for Penultimate (MB or TS) versus Ultimate (HMB or HTS) genes within the four standard, mutually exclusive gene categories. Uses the existing protein-coding reference, `SBGEandSSSdataForMBE.csv` (`rmf`), `Assoc_Data.txt` (`max(abs(Wald_Stat))` per mapped gene), and `sa_resolution_feature_membership.csv`; run the reference and SA-resolution wrangling first. A gene positive for either defining feature enters the corresponding stage once; genes positive for both stages remain in both. Missing feature flags do not count as positive. Metric-specific available-gene universes are retained rather than intersected. Dots show medians of 5,000 bootstrap means with 95% percentile intervals; grey lines show all-category arithmetic means without stage filtering. Exports `FigS7_rmf_*`, `FigS7_sa_effect_*`, and `FigS7_genome_wide_reference_means.csv` to `Data/Summary_tables/`. Each metric has four two-sided overlap-aware tests (10,000 permutations), with BH correction separately across its four tests. The former radar generator and outputs are archived in `side/archived_figS7_radar_context/`.

6. `figureS8_fly_human_conservation.R`
   Builds `Sup_figures/FigS8_fly_human_conservation.png` from the tracked compact outputs of the local conservation workflow in `fly_human_conservation/`. Supply the completed alignment and mapping directories explicitly, for example:

   ```bash
   Rscript Sup_plotting_and_analyses/figureS8_fly_human_conservation.R \
     --alignment-dir fly_human_conservation/alignment_<run-id> \
     --mapping-dir fly_human_conservation/mapping_<run-id>
   ```

   Exports category eligibility/missingness, bootstrap summaries, and pairwise permutation tests as `FigS8_*.csv` in `Data/Summary_tables/`. Preserve both explicit directory arguments so the plotted inputs and mapping attrition are identifiable; see `fly_human_conservation/README.md` for mapping, alignment, and provenance details.

7. `figureS9_conservation_adjusted_sa.R`
   Builds `Sup_figures/FigS9_conservation_adjusted_sa.png`, showing raw and
   conservation-adjusted means on exactly the same measured genes, using all
   orthologues with both SA and conservation measurements (one panel):

   ```bash
   Rscript Sup_plotting_and_analyses/figureS9_conservation_adjusted_sa.R \
     --alignment-dir fly_human_conservation/alignment_<run-id>
   ```

   The response is `max(abs(Wald_Stat))` over finite SNP statistics overlapping
   each protein-coding gene in `Data/inputs/Assoc_Data.txt`, matching Fig1.
   The signed input is a Wald z-statistic; absolute values are taken before
   the maximum, without squaring. Models are ordinary least squares, with category alone or
   category plus whole-alignment identity and both coverage measures
   (additive, common linear slopes; centered/scaled within each fit). Adjusted
   means predict at the pooled covariate means of the relevant analysis set.
   Background is excluded because its conservation is undefined, never zero.

   Two hypotheses are tested before and after adjustment: a global category
   effect across the three orthologue categories, then the disease-only versus
   other-orthologues effect after removing SSDOs. No factorial set of pairwise
   residual tests is performed. BH adjustment covers these four tests within
   the full measured SA/conservation universe. `conservation_sa_regression.R` implements the inference;
   `Rscript Sup_plotting_and_analyses/test_conservation_sa_regression.R` checks
   it against ordinary regression and explicit cluster resampling.

   Shared human partners and unequal residual variances make simple gene-label
   or residual shuffling inappropriate. Instead, each null-model residual
   vector is multiplied by a random +1/-1 sign shared within each human-gene
   cluster; the full model is refitted 10,000 times and a CR1 cluster-robust
   Wald statistic compared with the observed value (plus-one P). This is an
   approximate restricted wild-cluster bootstrap, not an exact label
   permutation test. It allows heteroskedasticity and dependence within
   partner clusters, assuming clusters are independent and the mean model is
   adequate. The 95% percentile intervals resample human-gene clusters
   (5,000 replicates), refitting the model at fixed target covariate means.
   Algorithm reference: Roodman et al. (2019), *The Stata Journal* 19:4-60,
   [doi:10.1177/1536867X19830877](https://doi.org/10.1177/1536867X19830877);
   the [wild-cluster bootstrap overview](https://s3alfisc.github.io/fwildclusterboot/articles/wild_bootstrap.html)
   also gives the restricted residual-resampling and CR1 covariance formulas.

   Exports `FigS9_*.csv`: all-gene inclusion/missingness audit, measured-gene
   fitted values/residuals, summaries, regression coefficients and robust
   standard errors, tests, covariate centers/scales and model diagnostics,
   plus input checksums. `FigS9_R_session.txt` records the R environment.
   A residual-only category comparison would not generally equal covariate
   adjustment: the category design must also be residualized against the
   covariates (Frisch–Waugh–Lovell). The joint models avoid that problem.
   These observational adjustments assess attenuation, not causal mediation;
   unmeasured confounding, nonlinearity, and gene-family structure beyond the
   retained human-partner clusters remain limitations. The existing per-gene
   maximum is also sensitive to the number of mapped SNPs; `NSNPs` is retained
   in the audit, but SNP coverage is not an added covariate in this
   conservation-only follow-up.

The SBGE-stratified component PNGs are no longer kept as tracked standalone outputs.

`FigS1` and `FigS4` are generated by the corresponding scripts in `Plotting_and_analyses/`.
