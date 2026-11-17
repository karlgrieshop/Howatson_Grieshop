# Supplementary analyses

Run scripts from the repository root. Figures are written to `Sup_figures/`;
statistical tables and figure data are written to `Data/Summary_tables/`.

| Script | Output or analysis |
|---|---|
| `sbge_metric_plots.R` | Expression-stratified metric summaries |
| `figureS2_sa_mag.R` | Fig. S2 |
| `figureS3_context.R` | Fig. S3 and candidate-context summaries |
| `figureS5_figureS6_metric_panels.R` | Figs. S5–S6 |
| `figureS7_fly_human_conservation.R` | Fig. S7 and conservation tests |
| `figureS8_conservation_adjusted_sa.R` | Fig. S8 and clustered regression tests |
| `figureS9_resolution_stages.R` | Fig. S9 and stage-specific tests |

Fig. S7 requires explicit mapping and alignment directories. Fig. S8 requires
an alignment directory and `Data/inputs/Assoc_Data.txt`. Both use the
conservation methods and definitions described in
`fly_human_conservation/README.md`.
