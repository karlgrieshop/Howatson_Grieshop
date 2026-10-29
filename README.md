# SSDOs_Howatson_Grieshop_Phil_Trans

This repo includes the code to support Howatson and Grieshop (2026) Phil Trans B, a study of *Drosophila melanogaster* genes whose human orthologues are implicated in various sex-specific genetic diseases (i.e. sex-specific disease orthologues, SSDOs).

> **All raw inputs available from:**  
> **[ACCESSION TO INSERT]**
>
> The repository intentionally keeps `Data/inputs/` **gitignored**. To rerun the workflow, download or clone the repository and copy the accession files into the matching local paths under `Data/inputs/`.

## Contents

- **Main figures:** `Figures/Fig1_SA_Magnitude.png`, `Figures/Fig2_rmf_tau.png`, `Figures/Fig3_rlns_tajD_DoS.png`, `Figures/Fig4_sa_res_radar.png`
- **Supplementary figures:** `Sup_figures/FigS1_SA_Candidate_Enrichment.png`, `Sup_figures/FigS2_sa_mag.png`, `Sup_figures/FigS3_sa_mag_context.png`, `Sup_figures/FigS4_SA_Magnitude_Brigos5.png`, `Sup_figures/FigS5_rmf_tau.png`, `Sup_figures/FigS6_rlns_tajD_DoS.png`, `Sup_figures/FigS7_resolution_stages.png`
- **Tables and data frames:** `Data/outputs/`, `Data/Summary_tables/`, and generated Word tables in `Tables/`
- **Code:** `Data_wrangling/`, `Plotting_and_analyses/`, and `Sup_plotting_and_analyses/`

## Repository layout

```text
Data/
  inputs/          local accession files copied in by the user; intentionally untracked
  outputs/         canonical compiled reference/intermediate tables used by active scripts
  Summary_tables/  kept figure-driving summaries and retained statistical tables
Data_wrangling/                active wrangling/build scripts
Plotting_and_analyses/         active main-figure and main-text analysis scripts
Sup_plotting_and_analyses/     active supplementary figure scripts
Tables/            table-generation scripts and local Word outputs
Figures/           final main-figure PNGs (Fig1–Fig4)
Sup_figures/       final supplementary PNGs (FigS1–FigS7)
```

## Input provenance

All files below should be supplied locally in `Data/inputs/` from the accession above. Fill in the DOI placeholders before release.

| Local file(s) | Original source | DOI / link |
| --- | --- | --- |
| `SBGEandSSSdataForMBE.csv` | Singh and Agrawal (2023) | [add DOI] |
| `FileS7.txt` | Fraïsse et al. (2019) | [add DOI] |
| `Assoc_Data.txt`, `Ruzicka_et_al_pbio.3000244.s008.csv` | Ruzicka et al. (2019) | [add DOI] |
| `Innocenti_and_Morrow_ANTAGONISTIC.csv` | Innocenti and Morrow (2010) | [add DOI] |
| `DE.summary.tsv` | Grieshop et al. (2025); 695 unique FBgn genes with `Sig == TRUE` | [add DOI] |
| `DifferentialGeneExpression.whole.bodies.tsv`, `Sex.Reversed.AI.test.tsv`, `Heads.FDR.test.for.AI.tsv`, `Gonads.FDR.test.for.AI.tsv`, `Genes.with.Sex.Reversed.AI.tsv`, `Genes.for.analysis.final.tsv` | Mishra et al. ASE study | [add DOI] |
| `Mishra.sig.SSS_filt25.txt`, `dimorphic.subset.list_wDsRedgenes.txt`, `all.exclude.txt` | Grieshop et al. (2025) / MaleLimitedEvo locally supplied SSS tested, final, and excluded-gene lists | [add DOI / note] |
| `all_paralogs_correct.csv` | Harper and Morrow (2025) | [add DOI] |
| `Drosophila_melanogaster.BDGP6.54.115.gtf` | Ensembl / FlyBase gene annotation release 115 | [add link] |
| `all_disease_orthologs.csv`, `all_orthologs.csv`, `ssdo_inputs/**` | This study’s curated SSDO / orthology source files | [add accession note] |

## Fastest route: rebuild figures from kept summaries

The repository keeps the summary tables needed for the **already-summarized panel figures**. If those tables are present, readers can jump straight to:

1. `Rscript Plotting_and_analyses/figure1_sa_magnitude.R`
2. `Rscript Plotting_and_analyses/figure2_figure3_metric_panels.R`
3. `Rscript Sup_plotting_and_analyses/figureS5_figureS6_metric_panels.R`

This shortcut works directly from the kept plot-data tables in `Data/Summary_tables/`.

**Still upstream-dependent:** `FigS1`, `FigS2`, `FigS3`, `Fig4`, `FigS4`, and `FigS7` are rebuilt from their active analysis scripts rather than from preassembled panel-summary tables alone.

## Full active workflow

Run these from the repository root.

1. `Rscript Data_wrangling/compile_ssdo_reference_data.R`
2. `Rscript Data_wrangling/build_mishra_whole_body_sbge_bins.R`
3. `Rscript Plotting_and_analyses/boot_permute_metrics.R`
4. `Rscript Plotting_and_analyses/sa_magnitude_plots.R`
5. `Rscript Plotting_and_analyses/figure1_sa_magnitude.R`
6. `Rscript Plotting_and_analyses/figure2_figure3_metric_panels.R`
7. `Rscript Sup_plotting_and_analyses/sbge_metric_plots.R`
8. `Rscript Sup_plotting_and_analyses/figureS2_sa_mag.R`
9. `Rscript Sup_plotting_and_analyses/figureS3_context.R`
10. `Rscript Sup_plotting_and_analyses/figureS5_figureS6_metric_panels.R`
11. `Rscript Data_wrangling/build_affx_dros_2_array_reference.R` *(only if `Data/outputs/Affx_Dros_2.0.csv` needs to be rebuilt)*
12. `Rscript Plotting_and_analyses/sa_candidate_enrichment.R`
13. `Rscript Data_wrangling/build_sa_resolution_inputs.R`
14. `Rscript Plotting_and_analyses/figure4_sa_resolution.R`
15. `Rscript Sup_plotting_and_analyses/figureS7_resolution_stages.R`
16. `Rscript Plotting_and_analyses/figureS4_figure1_with_brigos_5class.R`
17. `Rscript Tables/table1_sa_magnitude_permutation_docx.R`
18. `Rscript Tables/tableS1_candidate_enrichment_docx.R`
19. `Rscript Tables/tableS2_candidate_context_docx.R`

## What each active area does

- `Data_wrangling/compile_ssdo_reference_data.R` compiles the canonical SSDO, disease-orthologue, orthologue, and gene-reference tables from the accession inputs.
- `Data_wrangling/build_mishra_whole_body_sbge_bins.R` rebuilds the active five-bin Mishra whole-body SBGE classification using the MaleLimitedEvo thresholds.
- `Plotting_and_analyses/boot_permute_metrics.R` writes the kept four-group metric plot-data and permutation-test tables for Tajima’s *D*, *r*<sub>*MF*</sub>, *pi*<sub>*rN*</sub>, tissue specificity, and *DoS*, with BH correction applied within each metric across the full six pairwise group comparisons.
- `Plotting_and_analyses/sa_magnitude_plots.R` writes the Ruzicka SA-effect summaries used in Fig1 and the six pairwise permutation tests used in Table1, with BH correction across those six comparisons. The retained metric is the maximum SNP-level `Wald_Stat` per protein-coding gene (a Wald chi-squared statistic, not a conventional effect-size estimate).
- `Plotting_and_analyses/sa_candidate_enrichment.R` generates **FigS1** and **TableS1**, testing enrichment of Innocenti and Morrow, Ruzicka et al., and Brigos-Baril et al. (2026) candidate genes among SSDOs, disease orthologues, and orthologues. The first two tests use their study-matched assay universes; Brigos candidates are tested in the protein-coding-gene universe. TableS1 reports raw two-sided Fisher P-values and BH-adjusted P-values across all nine tests. FigS1 asterisks indicate BH-adjusted P < 0.05 across that same nine-test family.
- `Sup_plotting_and_analyses/figureS3_context.R` generates **FigS3**, comparing SSDOs with Innocenti and Morrow (2010), Grieshop et al. (2025), and Ruzicka et al. (2019) candidates using per-gene maximum Ruzicka Wald statistics. Of the 695 significant Grieshop DE genes, 530 have estimates in the current protein-coding/SNP-mapped universe (121 are outside the protein-coding reference; 44 more lack mapped estimates). Groups may overlap. Original Fig1 non-SSDO bootstrap means remain blue/green/grey reference lines. Gene-level values, bootstrap summaries, reference lines, and the Grieshop inclusion audit are exported to `Data/Summary_tables/FigS3_*.csv`.
- `Tables/tableS2_candidate_context_docx.R` exports the three SSDO-versus-candidate SA-effect comparisons and builds TableS2 using the existing overlap-aware permutation procedure, with BH correction across these three comparisons.
- `Sup_plotting_and_analyses/figureS2_sa_mag.R` generates **FigS2**, the SA-effect comparison stratified by whole-body sex-biased gene expression.
- `Data_wrangling/build_sa_resolution_inputs.R` assembles the study-matched SA-resolution feature tables for **Fig4** and the radar supplements, storing `TRUE`/`FALSE` only within each mechanism's source-study universe and `NA` outside it; it also exports the corresponding feature-universe summary. `SSS` now uses the locally supplied `Mishra.sig.SSS_filt25.txt` tested-gene frame together with the DsRed-restored final list in `dimorphic.subset.list_wDsRedgenes.txt`; `Dup` remains the only protein-coding-genome fallback denominator.
- `Plotting_and_analyses/figure4_sa_resolution.R` generates **Fig4** and the kept SA-resolution statistical summaries using those study-matched per-feature denominators.
- `Sup_plotting_and_analyses/figureS7_resolution_stages.R` generates **FigS7**: rMF (A) and SA effects (B) for Penultimate (MB or TS) versus Ultimate (HMB or HTS) genes in the four mutually exclusive gene categories. Grey lines are arithmetic means across all protein-coding genes with finite values for each metric, including genes outside either stage. The metric-specific universes are not restricted to genes measured for both metrics. Dots and whiskers show bootstrap medians of means and 95% percentile intervals (5,000 resamples). Each metric has four within-category, two-sided, overlap-aware permutation tests (10,000 permutations), with shared genes fixed in both groups and only exclusive memberships permuted. BH correction is applied separately within each metric's four tests. Gene values, inclusion audits, bootstrap summaries, reference means, and tests are exported as `Data/Summary_tables/FigS7_*.csv`. The former radar-overlay figure, generator, and multivariate test are archived in `side/archived_figS7_radar_context/`.
- `Plotting_and_analyses/figureS4_figure1_with_brigos_5class.R` generates the single-panel **FigS4** Brigos-Baril et al. (2026) sensitivity figure from the Ruzicka SA-effect metric.

## Notes

- The default gene universe for the kept active workflow is the **protein-coding FBgn gene set**.
- `Data/outputs/` contains active intermediates used by downstream scripts; `Data/Summary_tables/` contains retained figure summaries and kept statistical tables.
