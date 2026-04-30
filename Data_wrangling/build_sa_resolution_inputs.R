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
dimorphic_subset_path <- prefer_local_input(
  "dimorphic.subset.list.txt",
  male_limited_path("Analyses", "DS_analysis", "dimorphic.subset.list.txt")
)
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


# Tissue-specific expression is handled in three stages:
# 1. gate on high tau (top quartile genome-wide),
# 2. within that TSE subset, identify genes in the extreme male or female
#    tissue-PC tails,
# 3. define SS_TSE as asymmetric membership between the male and female gates.
# Tissue-specific expression is defined in three stages:
# 1. gate on high tau across the protein-coding genome,
# 2. within that TSE subset, detect genes in the male or female PC tails,
# 3. classify SS_TSE as asymmetric membership between the male and female gates.
load_tse_sets <- function(gene_universe_ids) {
  singh <- read.csv(input_path("SBGEandSSSdataForMBE.csv"), stringsAsFactors = FALSE)
  colnames(singh)[2] <- "FlyBaseID"

  tissue_metrics <- singh |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      tau = as.numeric(tao.FML.FA2),
      m.Comp.1 = as.numeric(`m.Comp.1`),
      m.Comp.2 = as.numeric(`m.Comp.2`),
      m.Comp.3 = as.numeric(`m.Comp.3`),
      f.Comp.1 = as.numeric(`f.Comp.1`),
      f.Comp.2 = as.numeric(`f.Comp.2`),
      f.Comp.3 = as.numeric(`f.Comp.3`)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), FlyBaseID %in% gene_universe_ids) |>
    distinct(FlyBaseID, .keep_all = TRUE)

  tau_cutoff <- unname(stats::quantile(tissue_metrics$tau, probs = 0.75, na.rm = TRUE))
  tse_ids <- tissue_metrics$FlyBaseID[!is.na(tissue_metrics$tau) & is.finite(tissue_metrics$tau) & tissue_metrics$tau >= tau_cutoff]
  tse_data <- tissue_metrics |>
    filter(FlyBaseID %in% tse_ids)

  # A gene is flagged if any male/female tissue PC places it in the extreme
  # tail. This makes the category sensitive to strongly unusual tissue-bias
  # structure on any of the leading tissue-specific axes.
  # A gene is flagged if any of the first three sex-specific tissue PCs places
  # it in the extreme tail. This keeps the definition sensitive to strongly
  # unusual tissue-expression structure on any leading axis.
  build_tail_flags <- function(data, prefix, lower_prob, upper_prob) {
    pcs <- paste0(prefix, c(".Comp.1", ".Comp.2", ".Comp.3"))
    thresholds <- bind_rows(lapply(pcs, function(pc_name) {
      values <- data[[pc_name]]
      values <- values[!is.na(values) & is.finite(values)]

      data.frame(
        Variable = pc_name,
        Lower = unname(stats::quantile(values, probs = lower_prob, na.rm = TRUE)),
        Upper = unname(stats::quantile(values, probs = upper_prob, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))

    flags <- Reduce(`|`, lapply(pcs, function(pc_name) {
      lower <- thresholds$Lower[thresholds$Variable == pc_name]
      upper <- thresholds$Upper[thresholds$Variable == pc_name]
      values <- data[[pc_name]]
      (!is.na(values)) & is.finite(values) & (values <= lower | values >= upper)
    }))

    list(flags = flags, thresholds = thresholds)
  }

  male_tail <- build_tail_flags(tse_data, prefix = "m", lower_prob = 0.025, upper_prob = 0.975)
  female_tail <- build_tail_flags(tse_data, prefix = "f", lower_prob = 0.025, upper_prob = 0.975)

  tse_data <- tse_data |>
    mutate(
      M_TSE = male_tail$flags,
      F_TSE = female_tail$flags,
      SS_TSE = xor(M_TSE, F_TSE)
    )

  tail_setting <- "2.5% tails (95% bulk retained)"
  if (sum(tse_data$M_TSE) == 0 || sum(tse_data$F_TSE) == 0 || sum(tse_data$SS_TSE) == 0) {
    male_tail <- build_tail_flags(tse_data, prefix = "m", lower_prob = 0.05, upper_prob = 0.95)
    female_tail <- build_tail_flags(tse_data, prefix = "f", lower_prob = 0.05, upper_prob = 0.95)

    tse_data <- tse_data |>
      mutate(
        M_TSE = male_tail$flags,
        F_TSE = female_tail$flags,
        SS_TSE = xor(M_TSE, F_TSE)
      )

    tail_setting <- "5% tails (90% bulk retained)"
  }

  threshold_summary <- bind_rows(
    data.frame(
      Step = "TSE gate",
      Variable = "tau",
      Lower = tau_cutoff,
      Upper = NA_real_,
      Rule = "top quartile of protein-coding genes",
      TailSetting = NA_character_,
      stringsAsFactors = FALSE
    ),
    male_tail$thresholds |>
      mutate(Step = "M_TSE gate", Rule = "outside lower/upper tail thresholds", TailSetting = tail_setting),
    female_tail$thresholds |>
      mutate(Step = "F_TSE gate", Rule = "outside lower/upper tail thresholds", TailSetting = tail_setting)
  ) |>
    select(Step, Variable, Lower, Upper, Rule, TailSetting)

  list(
    TSE = sort(unique(tse_ids)),
    M_TSE = sort(unique(tse_data$FlyBaseID[tse_data$M_TSE])),
    F_TSE = sort(unique(tse_data$FlyBaseID[tse_data$F_TSE])),
    SS_TSE = sort(unique(tse_data$FlyBaseID[tse_data$SS_TSE])),
    thresholds = threshold_summary,
    tail_setting = tail_setting
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


load_dr_set <- function() {
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
tse_sets <- load_tse_sets(gene_universe$FlyBaseID)
retained_ai_ids <- load_mishra_retained_ids()
ai_sets <- load_ai_sets(retained_ai_ids)
dr_set <- load_dr_set()
sSS_set <- load_sss_set()
dup_set <- load_dup_set()


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
  TSE = unique(tse_sets$TSE),
  M_TSE = unique(tse_sets$M_TSE),
  F_TSE = unique(tse_sets$F_TSE),
  SS_TSE = unique(tse_sets$SS_TSE),
  Dup = unique(dup_set),
  SSS = unique(sSS_set),
  ASE = unique(ai_sets$ASE),
  SS_ASE = unique(ai_sets$SS_ASE),
  TS_ASE = unique(ai_sets$TS_ASE),
  DR = unique(dr_set)
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
  Feature = c("HFB", "FB", "MB", "HMB", "TSE", "M_TSE", "F_TSE", "SS_TSE", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "DR"),
  Description = c(
    "Highly female-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules",
    "Female-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules (Moderately female-biased bin)",
    "Male-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules (Moderately male-biased bin)",
    "Highly male-biased whole-body SBGE genes from Mishra rebuilt with the MaleLimitedEvo x1 = 1, x2 = 5, y0 = 1, lfcSE rules",
    "Top-quartile tissue-specific genes based on Singh and Agrawal whole-body tau estimates",
    "TSE genes in the extreme male-PC tails of tissue specificity (m.Comp.1-3)",
    "TSE genes in the extreme female-PC tails of tissue specificity (f.Comp.1-3)",
    "Genes in M-TSE or F-TSE, but not both",
    "Genes with one or more paralogues in the Harper and Morrow duplication analysis",
    "Consistently sex-specific splicing gene set from Grieshop et al. (2025), derived using Mishra et al. (2024a) reference splicing profiles and filtered with three additional Mishra et al. (2024b) populations",
    "Genes with allelic imbalance in the filtered Mishra whole-body sex-reversed AI table (males or females significant)",
    "Genes with significant whole-body ASE in exactly one sex only in the filtered Mishra sex-reversed AI table",
    "Genes with significant ASE in exactly one of heads or gonads in the retained Mishra ASE tables",
    "Genes with dominance-reversed allelic imbalance from Genes.with.Sex.Reversed.AI.tsv"
  ),
  Source = c(
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/outputs/mishra_whole_body_sbge_bins.csv",
    "Data/inputs/SBGEandSSSdataForMBE.csv (tau top quartile)",
    "Data/inputs/SBGEandSSSdataForMBE.csv (male tissue PCs within TSE genes)",
    "Data/inputs/SBGEandSSSdataForMBE.csv (female tissue PCs within TSE genes)",
    "Derived from M_TSE and F_TSE",
    "Data/inputs/all_paralogs_correct.csv (local copy preferred)",
    "Data/inputs/dimorphic.subset.list.txt (local copy preferred)",
    "Data/inputs/Sex.Reversed.AI.test.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/Sex.Reversed.AI.test.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/{Heads,Gonads}.FDR.test.for.AI.tsv filtered to Genes.for.analysis.final.tsv (local copy preferred)",
    "Data/inputs/Genes.with.Sex.Reversed.AI.tsv"
  ),
  Status = c(rep("ready", 14)),
  GeneCount = unname(feature_counts[c("HFB", "FB", "MB", "HMB", "TSE", "M_TSE", "F_TSE", "SS_TSE", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "DR")]),
  Notes = c(
    rep(NA_character_, 4),
    sprintf("TSE is defined here as the top quartile of tau across the protein-coding gene universe (tau >= %.6f).", tse_sets$thresholds$Lower[tse_sets$thresholds$Step == "TSE gate"][1]),
    sprintf("Male tissue-PC tails were defined within TSE genes using %s.", tse_sets$tail_setting),
    sprintf("Female tissue-PC tails were defined within TSE genes using %s.", tse_sets$tail_setting),
    "SS_TSE is defined as genes in M_TSE or F_TSE, but not both.",
    "Dup is defined here as binary paralogue presence: genes that appear as FBgn_ID in the shared Harper and Morrow all_paralogs_correct.csv table. The current shared table already satisfies the Harper and Morrow paralogue criterion and contains only rows with DIOPT_score >= 2.",
    sprintf("Using dimorphic.subset.list.txt as the final consistently sex-specific splicing (SSS) set from Grieshop et al. (2025). Those 2036 source genes were defined from significant SSS genes in Mishra et al. (2024a), then retained only if male and female splicing profiles remained concordant with the Mishra reference profiles across three additional Mishra et al. (2024b) populations. Of those, %d were present in the current protein-coding analysis universe.", as.integer(feature_counts["SSS"])),
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
  tse_sets$thresholds,
  summary_path("SA_resolution_tse_thresholds.csv"),
  row.names = FALSE
)
