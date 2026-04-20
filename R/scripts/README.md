# Scripts

## Analysis Pipeline Order

Run scripts in this order for a full end-to-end analysis:

```
1. Final_Script.R                          → generates SBGE_disease_data.csv, mish_fra_SBGE.csv
2. data_analysis.R                         → generates all Summary Table CSVs
3. Candidate vs Ortholog and Disease genes.R → 4-way comparison plots
4. Plotting (Genomic Metrics, All).R       → standalone figure plots (reads Summary Table CSVs)
5. SA magnitude.R                          → sexually antagonistic SNP annotation (optional)
```

> All scripts use **absolute paths** hardcoded to the developer's local machine (`/Users/elliehowatson/...`). Update all `read.csv()`, `read.delim()`, and `setwd()` calls to match your local paths before running.

---

## Script Descriptions

### `Final_Script.R`

**Purpose**: Main pipeline — imports SBGE and population-genetics reference data, builds disease and background gene sets, and computes two-group (Disease vs Background) bootstrap summaries.

**Key steps**:
1. Load sex-specific disease gene lists (FlyBase disease-name searches + literature orthologue files)
2. Load SBGE data (Mishra 2024 whole-body), Fraisse 2018 (tissue-specificity / MK stats), Singh & Agrawal 2023 (Tajima's D, rMF, πN)
3. Define background (control) gene set = chromosomal genes minus disease genes
4. Merge disease and background genes with each reference dataset
5. Bootstrap 95% CIs and produce summary tables (saved as intermediate CSVs)

**Known fixed bugs**:
- `bind_rows(sexspecific_disease_data)` (was incorrectly `sexspecific_disease`; line 34)
- `select()+rename()` replaces `cbind()` for `SBG_whole_bodies` to preserve numeric types (line 70)
- `background_genes` now excludes disease genes (previously included all Fra genes; line 310)
- Chromosome files Xchr/Ychr now read with `header=FALSE` (previously `header=TRUE` dropped first gene; lines 193/197)
- Tissue-specificity background now explicitly filtered to `Tissue == "Whole Body"` with a `stopifnot()` assertion to prevent multi-tissue inflation (lines 377–390)

---

### `data_analysis.R`

**Purpose**: Extended analysis — multi-tissue SBGE, full summary-table generation for the 4-way comparison.

**Key steps**:
1. Same data loading as `Final_Script.R`
2. Builds multi-tissue SBGE dataset (`SBG_data` = gonads + heads + whole bodies)
3. Generates bootstrap summary tables for all metrics
4. Tissue-specificity analysis uses **whole-body SBGE only** (not multi-tissue) for background

**Known fixed bugs**:
- Xchr/Ychr `header=TRUE` → `header=FALSE`
- `disease` column removed from `distinct()` in tissue-specificity, DoS, rlnd comparisons (previously inflated multi-disease gene counts)
- `background_genes` now excludes disease genes
- `control_mish_SBGE_wb` filtered to `Tissue == "Whole Body"` + `stopifnot()` assertion

---

### `Candidate vs Ortholog and Disease genes.R`

**Purpose**: Assembles the four-way comparison (SSDOs / All Disease Orthologues / All Orthologues / Background), runs permutation tests, and produces the main comparison plots.

**Key steps**:
1. Loads pre-merged disease gene CSVs (from `Final_Script.R`)
2. Loads all-disease and all-ortholog gene lists, builds Fra and Singh merged datasets per group
3. Runs `TwoPerm()` and `TwoBoot()` permutation/bootstrap tests for each metric × SBGE bin
4. Produces multi-panel dot-whisker plots

**Known fixed bugs**:
- `GeneType = "Orthologs"` → `"Ortholog"` (was inconsistent between Fra and Singh datasets, breaking 4-way factor levels)

**Statistical notes**:
- Permutation tests are two-group label-permutation tests (disease vs. non-disease within the same dataset)
- Bootstrap uses stratified resampling within sex-bias category
- No FDR correction is applied across the many pairwise tests; interpret per-category p-values accordingly

---

### `Plotting (Genomic Metrics, All).R`

**Purpose**: Standalone plotting script. Reads pre-computed Summary Table CSVs and generates publication-quality figures.

**Figures produced**:
- `Fig1A`: Relative level of nonsynonymous diversity (rlnd / π̃ᵣₙ)
- `Fig1B`: Direction of Selection (DoS)
- `Fig2A`: Male-to-female fitness ratio (rMF)
- `Fig2B`: Tajima's *D*
- `Fig2C`: Tissue-specificity index

**Plot improvements applied**:
- Colorblind-safe palette (Wong 2011): SSDOs = `#D55E00`, Disease = `#0072B2`, Ortholog = `#009E73`, Background = `#999999`
- Shared `pub_theme()` helper for consistent publication formatting
- Removed duplicate `geom_point` + `geom_errorbar` layers (each plot previously drew error bars twice)

---

### `SA magnitude.R`

**Purpose**: Maps sexually antagonistic (SA) SNPs from a CMH test to genes via BioMart chromosomal-region queries.

**Dependencies**: `biomaRt`, `readxl`, `data.table`

**Output**: `gene_hits_snps.csv`, `gene_SA_magnitude.csv`
