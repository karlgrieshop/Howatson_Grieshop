# Assemble the binary gene-by-mechanism table used by the core SA-resolution
# analyses. Each row is one protein-coding FlyBase gene; each mechanism column
# is a logical TRUE/FALSE indicator showing whether that gene qualifies for the
# corresponding putative resolution mechanism.

# Assemble the binary gene-by-mechanism matrix used by the core SA-resolution
# analyses. Each row is one protein-coding FlyBase gene and each mechanism
# column is a logical TRUE/FALSE indicator for a putative route of conflict
# resolution.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Data_wrangling")) &&
      file.exists(file.path(current, "README.md"))
    ) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate the project root. Run this script from inside the repository root or one of its subdirectories.",
        call. = FALSE
      )
    }

    current <- parent
  }
}


project_root <- find_project_root()
input_path <- function(...) file.path(project_root, "Data", "inputs", ...)
output_path <- function(...) file.path(project_root, "Data", "outputs", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
mishra_ase_path <- function(...) file.path(project_root, "..", "Mishra_et_al_Drosophila_ASE", ...)
shared_parent_path <- function(...) file.path(project_root, "..", ...)
male_limited_path <- function(...) file.path(project_root, "..", "MaleLimitedEvo", ...)
harper_dup_path <- function(...) {
  file.path(
    dirname(dirname(dirname(project_root))),
    "BIO UEA",
    "Grieshop Lab",
    "Grieshop Lab shared",
    "Data",
    "Harper and Morrow 2025",
    ...
  )
}


suppressPackageStartupMessages({
  library(dplyr)
})


# Prefer local copies in `Data/inputs/` when present so the manuscript workflow
# can be rerun from the repo surface without depending on sibling directories.
# Prefer local untracked copies in `Data/inputs/` so the publication workflow
# can be rerun from the repo surface without depending on sibling directories.
prefer_local_input <- function(local_name, fallback_path) {
  local_path <- input_path(local_name)
  if (file.exists(local_path)) local_path else fallback_path
}


genes_for_analysis_path <- prefer_local_input(
  "Genes.for.analysis.final.tsv",
  shared_parent_path("Genes.for.analysis.final.tsv")
)
sex_reversed_ai_path <- prefer_local_input(
  "Sex.Reversed.AI.test.tsv",
  mishra_ase_path("Data_for_R_scripts", "Sex.Reversed.AI", "Sex.Reversed.AI.test.tsv")
)
heads_fdr_path <- prefer_local_input(
  "Heads.FDR.test.for.AI.tsv",
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Heads.FDR.test.for.AI.tsv")
)
gonads_fdr_path <- prefer_local_input(
  "Gonads.FDR.test.for.AI.tsv",
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Gonads.FDR.test.for.AI.tsv")
)
dimorphic_subset_path <- input_path("dimorphic.subset.list.txt")
duplication_table_path <- prefer_local_input(
  "all_paralogs_correct.csv",
  harper_dup_path("Output", "all_paralogs_correct.csv")
)


required_files <- c(
  output_path("dmel_gene_reference.csv"),
  output_path("sexspecific_disease_ids.csv"),
  output_path("mishra_whole_body_sbge_bins.csv"),
  input_path("SBGEandSSSdataForMBE.csv"),
  genes_for_analysis_path,
  sex_reversed_ai_path,
  heads_fdr_path,
  gonads_fdr_path,
  input_path("Genes.with.Sex.Reversed.AI.tsv"),
  dimorphic_subset_path,
  duplication_table_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required files for SA-resolution input assembly:\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


dir.create(output_path(), recursive = TRUE, showWarnings = FALSE)
dir.create(summary_path(), recursive = TRUE, showWarnings = FALSE)


# The analysis universe is explicitly the current protein-coding FlyBase gene
# set parsed from the GTF-derived reference table, not every valid-looking FBgn
# that appears in upstream source files.
# The analysis universe is explicitly the current protein-coding FlyBase gene
# set parsed from the GTF-derived reference table.
load_gene_universe <- function() {
  read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE) |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      Fly.Symbol = as.character(Fly.Symbol),
      gene_biotype = as.character(gene_biotype),
      chromosome = as.character(chromosome)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding") |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    arrange(FlyBaseID)
}


load_ssdo_ids <- function() {
  read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
    pull(FlyBaseID) |>
    as.character() |>
    unique()
}


# Rebuild the Mishra whole-body SBGE bins as four binary mechanism sets.
# Re-express the Mishra SBGE bins as four binary mechanism sets used downstream.
load_sbge_sets <- function() {
  sbge <- read.csv(output_path("mishra_whole_body_sbge_bins.csv"), stringsAsFactors = FALSE) |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      Category = as.character(Sex.Bias.Category)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), !is.na(Category)) |>
    distinct(FlyBaseID, .keep_all = TRUE)

  list(
    HFB = sbge$FlyBaseID[sbge$Category == "Highly female-biased"],
    FB = sbge$FlyBaseID[sbge$Category == "Moderately female-biased"],
    MB = sbge$FlyBaseID[sbge$Category == "Moderately male-biased"],
    HMB = sbge$FlyBaseID[sbge$Category == "Highly male-biased"]
  )
}


# The current Fig4 workflow no longer uses tissue-PC tail subsets. Instead it
# partitions the protein-coding gene universe into tau quintiles and retains
# four reader-facing categories: highest tissue specificity (HTS), the next-most
# tissue-specific quintile (TS), the next-most tissue-general quintile (TG), and
# highest tissue generality (HTG). The middle quintile is recorded in the
# threshold table for transparency but not carried into the radar.
load_tau_quintile_sets <- function(gene_universe_ids) {
  singh <- read.csv(input_path("SBGEandSSSdataForMBE.csv"), stringsAsFactors = FALSE)
  colnames(singh)[2] <- "FlyBaseID"

  tissue_metrics <- singh |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      tau = as.numeric(tao.FML.FA2)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), FlyBaseID %in% gene_universe_ids, !is.na(tau), is.finite(tau)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    mutate(TauQuintile = dplyr::ntile(tau, 5))

  quintile_ranges <- tissue_metrics |>
    group_by(TauQuintile) |>
    summarise(
      Lower = min(tau, na.rm = TRUE),
      Upper = max(tau, na.rm = TRUE),
      GeneCount = dplyr::n(),
      .groups = "drop"
    )

  threshold_summary <- data.frame(
    Step = c("HTG gate", "TG gate", "NA middle quintile", "TS gate", "HTS gate"),
    Variable = "tau",
    Lower = quintile_ranges$Lower[match(c(1, 2, 3, 4, 5), quintile_ranges$TauQuintile)],
    Upper = quintile_ranges$Upper[match(c(1, 2, 3, 4, 5), quintile_ranges$TauQuintile)],
    Rule = c(
      "lowest tau quintile (bottom 20 percent)",
      "second-lowest tau quintile",
      "middle tau quintile (not plotted)",
      "second-highest tau quintile",
      "highest tau quintile (top 20 percent)"
    ),
    TailSetting = NA_character_,
    GeneCount = quintile_ranges$GeneCount[match(c(1, 2, 3, 4, 5), quintile_ranges$TauQuintile)],
    stringsAsFactors = FALSE
  ) |>
    select(Step, Variable, Lower, Upper, Rule, TailSetting)

  list(
    HTG = sort(unique(tissue_metrics$FlyBaseID[tissue_metrics$TauQuintile == 1])),
    TG = sort(unique(tissue_metrics$FlyBaseID[tissue_metrics$TauQuintile == 2])),
    TS = sort(unique(tissue_metrics$FlyBaseID[tissue_metrics$TauQuintile == 4])),
    HTS = sort(unique(tissue_metrics$FlyBaseID[tissue_metrics$TauQuintile == 5])),
    thresholds = threshold_summary
  )
}


load_mishra_retained_ids <- function() {
  read.delim(genes_for_analysis_path, stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(geneID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct() |>
    pull(FlyBaseID)
}


# These ASE-derived mechanism classes all start from the retained Mishra gene
# set so that the manuscript uses the same filtering frame as the source study.
read_mishra_sex_reversed_ai <- function(retained_ids) {
  read.delim(
    sex_reversed_ai_path,
    stringsAsFactors = FALSE
  ) |>
    transmute(
      FlyBaseID = as.character(geneID),
      MaleAI = !is.na(p.val.males) & as.numeric(p.val.males) < 0.05,
      FemaleAI = !is.na(p.val.females) & as.numeric(p.val.females) < 0.05
    ) |>
    filter(
      !is.na(FlyBaseID),
      grepl("^FBgn", FlyBaseID),
      FlyBaseID %in% retained_ids
    ) |>
    distinct(FlyBaseID, .keep_all = TRUE)
}


read_mishra_sex_ai_fdr <- function(file_name, retained_ids) {
  read.delim(
    file_name,
    stringsAsFactors = FALSE
  ) |>
    transmute(
      FlyBaseID = as.character(geneID),
      HasAI = !is.na(FDR.intercept) & as.numeric(FDR.intercept) < 0.05,
      HasSexDependentAI = !is.na(FDR.sex.term) & as.numeric(FDR.sex.term) < 0.05
    ) |>
    filter(
      !is.na(FlyBaseID),
      grepl("^FBgn", FlyBaseID),
      FlyBaseID %in% retained_ids
    ) |>
    distinct(FlyBaseID, .keep_all = TRUE)
}


# ASE-derived mechanism classes are restricted to the retained Mishra gene set
# so the manuscript inherits the source-study filtering frame.
load_ai_sets <- function(retained_ids) {
  whole_body_ai <- read_mishra_sex_reversed_ai(retained_ids)
  heads <- read_mishra_sex_ai_fdr(heads_fdr_path, retained_ids)
  gonads <- read_mishra_sex_ai_fdr(gonads_fdr_path, retained_ids)

  head_gonad_ai <- full_join(
    heads |>
      select(FlyBaseID, HeadsAI = HasAI),
    gonads |>
      select(FlyBaseID, GonadsAI = HasAI),
    by = "FlyBaseID"
  ) |>
    mutate(
      HeadsAI = ifelse(is.na(HeadsAI), FALSE, HeadsAI),
      GonadsAI = ifelse(is.na(GonadsAI), FALSE, GonadsAI)
    )

  list(
    ASE = sort(unique(whole_body_ai$FlyBaseID[whole_body_ai$MaleAI | whole_body_ai$FemaleAI])),
    SS_ASE = sort(unique(whole_body_ai$FlyBaseID[xor(whole_body_ai$MaleAI, whole_body_ai$FemaleAI)])),
    TS_ASE = sort(unique(head_gonad_ai$FlyBaseID[xor(head_gonad_ai$HeadsAI, head_gonad_ai$GonadsAI)]))
  )
}


load_ssdr_set <- function() {
  read.delim(input_path("Genes.with.Sex.Reversed.AI.tsv"), stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(geneID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


load_sss_set <- function() {
  read.delim(
    dimorphic_subset_path,
    header = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transmute(FlyBaseID = as.character(V1)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


load_dup_set <- function() {
  read.csv(duplication_table_path, stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(FBgn_ID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


gene_universe <- load_gene_universe()
ssdo_ids <- load_ssdo_ids()
sbge_sets <- load_sbge_sets()
tau_sets <- load_tau_quintile_sets(gene_universe$FlyBaseID)
retained_ai_ids <- load_mishra_retained_ids()
ai_sets <- load_ai_sets(retained_ai_ids)
ssdr_set <- load_ssdr_set()
sSS_set <- load_sss_set()
dup_set <- load_dup_set()
sss_source_count <- length(unique(sSS_set))


# The final feature lists can be broader than the protein-coding universe in
# their raw source files, so the actual membership table below is treated as the
# source of truth for exported counts and feature-set IDs.
# Raw source lists can extend beyond the final protein-coding universe, so the
# exported counts below are derived from `feature_membership`, not the source
# vectors themselves.
available_feature_sets <- list(
  HFB = unique(sbge_sets$HFB),
  FB = unique(sbge_sets$FB),
  MB = unique(sbge_sets$MB),
  HMB = unique(sbge_sets$HMB),
  HTS = unique(tau_sets$HTS),
  TS = unique(tau_sets$TS),
  TG = unique(tau_sets$TG),
  HTG = unique(tau_sets$HTG),
  Dup = unique(dup_set),
  SSS = unique(sSS_set),
  ASE = unique(ai_sets$ASE),
  SS_ASE = unique(ai_sets$SS_ASE),
  TS_ASE = unique(ai_sets$TS_ASE),
  SSDR = unique(ssdr_set)
)


# ---- Main execution block -------------------------------------------------
# Build the analysis matrix, export it in both wide and long formats, then
# write a reader-facing notes table that records provenance and thresholds.

# ---- Main assembly workflow -----------------------------------------------

feature_membership <- gene_universe |>
  mutate(IsSSDO = FlyBaseID %in% ssdo_ids)

for (feature_name in names(available_feature_sets)) {
  feature_membership[[feature_name]] <- feature_membership$FlyBaseID %in% available_feature_sets[[feature_name]]
}


feature_sets_long <- bind_rows(lapply(names(available_feature_sets), function(feature_name) {
  data.frame(
    Feature = feature_name,
    FlyBaseID = sort(unique(feature_membership$FlyBaseID[feature_membership[[feature_name]]])),
    stringsAsFactors = FALSE
  )
}))

feature_counts <- vapply(names(available_feature_sets), function(feature_name) {
  sum(feature_membership[[feature_name]], na.rm = TRUE)
}, numeric(1))


feature_notes <- data.frame(
  Feature = c("HFB", "FB", "MB", "HMB", "HTS", "TS", "TG", "HTG", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "SSDR"),
  Description = c(
    "Highly female-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules",
    "Female-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules (Moderately female-biased bin)",
    "Male-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules (Moderately male-biased bin)",
    "Highly male-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules",
    "Highest-tau quintile from Singh and Agrawal whole-body tau estimates (highest tissue specificity)",
    "Second-highest tau quintile from Singh and Agrawal whole-body tau estimates",
    "Second-lowest tau quintile from Singh and Agrawal whole-body tau estimates",
    "Lowest-tau quintile from Singh and Agrawal whole-body tau estimates (highest tissue generality)",
    "Genes with one or more paralogues in the Harper and Morrow duplication analysis",
    "Consistently sex-specific splicing gene set from Grieshop et al. (2025), derived using Mishra et al. (2024a) reference splicing profiles and filtered with three additional Mishra et al. (2024b) populations",
    "Genes with allelic imbalance in the filtered Mishra whole-body sex-reversed AI table (males or females significant)",
    "Genes with significant whole-body ASE in exactly one sex only in the filtered Mishra sex-reversed AI table",
    "Genes with significant ASE in exactly one of heads or gonads in the retained Mishra ASE tables",
    "Genes with sex-specific dominance-reversed allelic imbalance from Genes.with.Sex.Reversed.AI.tsv"
  ),
  Source = c(
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/inputs/SBGEandSSSdataForMBE.csv (tau highest quintile)",
    "Data/inputs/SBGEandSSSdataForMBE.csv (tau second-highest quintile)",
    "Data/inputs/SBGEandSSSdataForMBE.csv (tau second-lowest quintile)",
    "Data/inputs/SBGEandSSSdataForMBE.csv (tau lowest quintile)",
    "Data/inputs/all_paralogs_correct.csv (local copy preferred)",
    "Data/inputs/dimorphic.subset.list.txt",
    "Data/inputs/Sex.Reversed.AI.test.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/Sex.Reversed.AI.test.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/{Heads,Gonads}.FDR.test.for.AI.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/Genes.with.Sex.Reversed.AI.tsv"
  ),
  Status = c(rep("ready", 14)),
  GeneCount = unname(feature_counts[c("HFB", "FB", "MB", "HMB", "HTS", "TS", "TG", "HTG", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "SSDR")]),
  Notes = c(
    rep(NA_character_, 4),
    sprintf("HTS is the highest tau quintile, with tau in [%.6f, %.6f].", tau_sets$thresholds$Lower[tau_sets$thresholds$Step == "HTS gate"][1], tau_sets$thresholds$Upper[tau_sets$thresholds$Step == "HTS gate"][1]),
    sprintf("TS is the second-highest tau quintile, with tau in [%.6f, %.6f].", tau_sets$thresholds$Lower[tau_sets$thresholds$Step == "TS gate"][1], tau_sets$thresholds$Upper[tau_sets$thresholds$Step == "TS gate"][1]),
    sprintf("TG is the second-lowest tau quintile, with tau in [%.6f, %.6f].", tau_sets$thresholds$Lower[tau_sets$thresholds$Step == "TG gate"][1], tau_sets$thresholds$Upper[tau_sets$thresholds$Step == "TG gate"][1]),
    sprintf("HTG is the lowest tau quintile, with tau in [%.6f, %.6f].", tau_sets$thresholds$Lower[tau_sets$thresholds$Step == "HTG gate"][1], tau_sets$thresholds$Upper[tau_sets$thresholds$Step == "HTG gate"][1]),
    "Dup is defined here as binary paralogue presence: genes that appear as FBgn_ID in the shared Harper and Morrow all_paralogs_correct.csv table. The current shared table already satisfies the Harper and Morrow paralogue criterion and contains only rows with DIOPT_score >= 2.",
    sprintf("Using Data/inputs/dimorphic.subset.list.txt as the final consistently sex-specific splicing (SSS) set from Grieshop et al. (2025). This local file contains the %d-gene final SSS list, including genes restored after the earlier DsRed-proximity exclusion was removed. Of those, %d were present in the current protein-coding analysis universe.", sss_source_count, as.integer(feature_counts["SSS"])),
    rep(NA_character_, 4)
  ),
  stringsAsFactors = FALSE
)


write.csv(
  feature_membership,
  output_path("sa_resolution_feature_membership.csv"),
  row.names = FALSE
)

write.csv(
  feature_sets_long,
  output_path("sa_resolution_feature_sets.csv"),
  row.names = FALSE
)

write.csv(
  feature_notes,
  summary_path("sa_resolution_feature_notes.csv"),
  row.names = FALSE
)

write.csv(
  tau_sets$thresholds,
  summary_path("SA_resolution_tau_quantiles.csv"),
  row.names = FALSE
)
