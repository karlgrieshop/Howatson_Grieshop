# SSDOs_Howatson_Grieshop_Phil_Trans

This repo includes the code to support Howatson and Grieshop (2026) Phil Trans B, a study of *Drosophila melanogaster* genes whose human orthologues are implicated in various sex-specific genetic diseases (i.e. sex-specific disease orthologues, SSDOs).

> **All raw inputs available from:**  
> **[ACCESSION TO INSERT]**
>
> The repository intentionally keeps `Data/inputs/` **gitignored**. To rerun the workflow, download or clone the repository and copy the accession files into the matching local paths under `Data/inputs/`.

## Contents

- **Main figures:** `Figures/Fig1_SA_Probability_and_Magnitude.png`, `Figures/Fig2_rmf_tau.png`, `Figures/Fig3_rlns_tajD_DoS.png`, `Figures/Fig4_sa_res_radar.png`
- **Supplementary figures:** `Sup_figures/FigS1_SA_Candidate_Enrichment.png`, `Sup_figures/FigS2_sa_prob_mag.png`, `Sup_figures/FigS3_rmf_tau.png`, `Sup_figures/FigS4_rlns_tajD_DoS.png`, `Sup_figures/FigS5_radar_context.png`, `Sup_figures/FigS6_SSDO_candidate_overlay.png`
- **Tables and data frames:** `Data/outputs/` and `Data/Summary_tables/`
- **Code:** `Data_wrangling/`, `Plotting/`, `Sup_plotting/`, and `Analyses/`

## Repository layout

```text
Data/
  inputs/          local accession files copied in by the user; intentionally untracked
  outputs/         canonical compiled reference/intermediate tables used by active scripts
  Summary_tables/  kept figure-driving summaries and retained statistical tables
Data_wrangling/    active wrangling/build scripts
Analyses/          active analysis scripts for FigS1 and Fig4
Plotting/          active main-figure plotting scripts
Sup_plotting/      active supplementary plotting scripts
Figures/           final main-figure PNGs (Fig1–Fig4)
Sup_figures/       final supplementary PNGs (FigS1–FigS6)
```

## Input provenance

All files below should be supplied locally in `Data/inputs/` from the accession above. Fill in the DOI placeholders before release.

| Local file(s) | Original source | DOI / link |
| --- | --- | --- |
| `SBGEandSSSdataForMBE.csv` | Singh and Agrawal (2023) | [add DOI] |
| `FileS7.txt` | Fraïsse et al. (2019) | [add DOI] |
| `Assoc_Data.txt`, `Ruzicka_et_al_pbio.3000244.s008.csv` | Ruzicka et al. (2019) | [add DOI] |
| `true_cmh_results.txt`, `gene_hits_snps_90.csv` | Neal, Truter, *et al.* (unpublished) | [add accession / note] |
| `Innocenti_and_Morrow_ANTAGONISTIC.csv` | Innocenti and Morrow (2010) | [add DOI] |
| `DifferentialGeneExpression.whole.bodies.tsv`, `Sex.Reversed.AI.test.tsv`, `Heads.FDR.test.for.AI.tsv`, `Gonads.FDR.test.for.AI.tsv`, `Genes.with.Sex.Reversed.AI.tsv`, `Genes.for.analysis.final.tsv` | Mishra et al. ASE study | [add DOI] |
| `dimorphic.subset.list.txt` | Grieshop et al. (2025) final SSS list, with the DsRed-restored genes included | [add DOI] |
| `all_paralogs_correct.csv` | Harper and Morrow (2025) | [add DOI] |
| `Drosophila_melanogaster.BDGP6.54.115.gtf` | Ensembl / FlyBase gene annotation release 115 | [add link] |
| `all_disease_orthologs.csv`, `all_orthologs.csv`, `ssdo_inputs/**` | This study’s curated SSDO / orthology source files | [add accession note] |

## Fastest route: rebuild figures from kept summaries

The repository keeps the summary tables needed for the **already-summarized panel figures**. If those tables are present, readers can jump straight to:

1. `Rscript Plotting/figure1_sa_probability_magnitude.R`
2. `Rscript Plotting/figure2_figure3_metric_panels.R`
3. `Rscript Sup_plotting/figureS3_figureS4_metric_panels.R`

This shortcut works directly from the kept plot-data tables in `Data/Summary_tables/`.

**Still upstream-dependent:** `FigS1`, `FigS2`, `Fig4`, `FigS5`, and `FigS6` are rebuilt from their active analysis scripts rather than from preassembled panel-summary tables alone.

## Full active workflow

Run these from the repository root.

1. `Rscript Data_wrangling/compile_ssdo_reference_data.R`
2. `Rscript Data_wrangling/build_mishra_whole_body_sbge_bins.R`
3. `Rscript Plotting/plotting.R`
4. `Rscript Plotting/sa_probability_plots.R`
5. `Rscript Plotting/sa_magnitude_plots.R`
6. `Rscript Plotting/figure1_sa_probability_magnitude.R`
7. `Rscript Plotting/figure2_figure3_metric_panels.R`
8. `Rscript Sup_plotting/sbge_metric_plots.R`
9. `Rscript Sup_plotting/figureS2_sa_prob_mag.R`
10. `Rscript Sup_plotting/figureS3_figureS4_metric_panels.R`
11. `Rscript Data_wrangling/build_affx_dros_2_array_reference.R` *(only if `Data/outputs/Affx_Dros_2.0.csv` needs to be rebuilt)*
12. `Rscript Analyses/sa_candidate_enrichment.R`
13. `Rscript Data_wrangling/build_sa_resolution_inputs.R`
14. `Rscript Analyses/core_analyses.R`
15. `Rscript Sup_plotting/figureS5_radar_context.R`

## What each active area does

- `Data_wrangling/compile_ssdo_reference_data.R` compiles the canonical SSDO, disease-orthologue, orthologue, and gene-reference tables from the accession inputs.
- `Data_wrangling/build_mishra_whole_body_sbge_bins.R` rebuilds the active five-bin Mishra whole-body SBGE classification using the MaleLimitedEvo thresholds.
- `Plotting/plotting.R` writes the kept four-group metric plot-data and permutation-test tables for Tajima’s *D*, *r*<sub>*MF*</sub>, *pi*<sub>*rN*</sub>, tissue specificity, and *DoS*, with BH correction applied within each metric across the full six pairwise group comparisons.
- `Plotting/sa_probability_plots.R` and `Plotting/sa_magnitude_plots.R` write the kept study-level SA summary tables used in Fig1.
- `Analyses/sa_candidate_enrichment.R` generates **FigS1** and **TableS1**.
- `Data_wrangling/build_sa_resolution_inputs.R` assembles the feature-membership table for **Fig4**, including whole-body SBGE classes, four non-middle tau quintiles (`HTS`, `TS`, `TG`, `HTG`), Dup, SSS, and Mishra-derived ASE classes.
- `Analyses/core_analyses.R` generates **Fig4** and the kept SA-resolution statistical summaries.
- `Sup_plotting/figureS5_radar_context.R` generates **FigS5**, a four-panel SA-resolution radar-context figure with panel A reprinting the Fig4 SSDO comparison and panels B-D replacing the orange focal group with the Innocenti, Ruzicka, and Neal candidate SA-gene sets while retaining the same split background categories; it also generates **FigS6**, a single-panel overlay comparing SSDOs directly with the three study-specific SA-candidate sets.

## Notes

- The default gene universe for the kept active workflow is the **protein-coding FBgn gene set**.
- `Data/outputs/` contains active intermediates used by downstream scripts; `Data/Summary_tables/` contains retained figure summaries and kept statistical tables.
