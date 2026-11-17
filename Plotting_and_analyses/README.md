# Main-text analyses

Run scripts from the repository root. Statistical summaries are written to
`Data/Summary_tables/`; figures are written to `Figures/`.

| Script | Output |
|---|---|
| `boot_permute_metrics.R` | Bootstrap summaries and permutation tests for Figs. 2–3 |
| `sa_magnitude_plots.R` | SA-effect summaries and tests for Fig. 1 and Table 1 |
| `figure1_sa_magnitude.R` | Fig. 1 |
| `figure2_figure3_metric_panels.R` | Figs. 2–3 |
| `sa_candidate_enrichment.R` | Fig. S1 and Table S1 |
| `figure4_sa_resolution.R` | Fig. 4 and SA-resolution tests |
| `figureS4_figure1_with_brigos_5class.R` | Fig. S4 |

The SA-effect metric is the per-gene maximum of absolute SNP-level
`Wald_Stat` values overlapping protein-coding genes. Absolute values are taken
before the maximum; values are not squared. Table-generation scripts are in
`Tables/`.
