# Data wrangling

Run scripts from the repository root. They write compiled inputs to
`Data/outputs/`.

| Script | Output or purpose |
|---|---|
| `compile_ssdo_reference_data.R` | Canonical gene categories and FlyBase gene/transcript references |
| `build_mishra_whole_body_sbge_bins.R` | Whole-body sex-biased expression bins |
| `build_affx_dros_2_array_reference.R` | Affymetrix probe-to-gene reference for the matched assay universe |
| `build_sa_resolution_inputs.R` | Study-matched feature sets and universes for Fig. 4 and Fig. S9 |

The feature-membership table records `TRUE`/`FALSE` within each source-study
universe and `NA` outside it. All source files required for these steps should
be supplied in `Data/inputs/`.
