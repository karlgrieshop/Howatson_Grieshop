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
        "Could not locate the project root. Run this script from inside the Final_Year_Project repository.",
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


required_files <- c(
  output_path("dmel_gene_reference.csv"),
  output_path("sexspecific_disease_ids.csv"),
  input_path("Sex.Biased.Genes.Whole.csv"),
  shared_parent_path("Genes.for.analysis.final.tsv"),
  mishra_ase_path("Data_for_R_scripts", "Sex.Reversed.AI", "Sex.Reversed.AI.test.tsv"),
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Whole.bodies.FDR.test.for.AI.tsv"),
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Reciprocal.whole.bodies.FDR.test.for.AI.tsv"),
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Heads.FDR.test.for.AI.tsv"),
  mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", "Gonads.FDR.test.for.AI.tsv"),
  input_path("Genes.with.Sex.Reversed.AI.tsv"),
  male_limited_path("Analyses", "DS_analysis", "dimorphic.subset.list.txt"),
  harper_dup_path("Output", "all_paralogs_correct.csv")
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


load_sbge_sets <- function() {
  sbge <- read.csv(input_path("Sex.Biased.Genes.Whole.csv"), stringsAsFactors = FALSE)
  colnames(sbge)[1] <- "FlyBaseID"

  sbge <- sbge |>
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


load_mishra_retained_ids <- function() {
  read.delim(shared_parent_path("Genes.for.analysis.final.tsv"), stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(geneID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct() |>
    pull(FlyBaseID)
}


read_mishra_sex_reversed_ai <- function(retained_ids) {
  read.delim(
    mishra_ase_path("Data_for_R_scripts", "Sex.Reversed.AI", "Sex.Reversed.AI.test.tsv"),
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
    mishra_ase_path("Data_for_R_scripts", "GLM.outputs", "FDR.Sex.Effects.in.AI", file_name),
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


load_ai_sets <- function(retained_ids) {
  whole_body_ai <- read_mishra_sex_reversed_ai(retained_ids)
  heads <- read_mishra_sex_ai_fdr("Heads.FDR.test.for.AI.tsv", retained_ids)
  gonads <- read_mishra_sex_ai_fdr("Gonads.FDR.test.for.AI.tsv", retained_ids)

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
    male_limited_path("Analyses", "DS_analysis", "dimorphic.subset.list.txt"),
    header = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transmute(FlyBaseID = as.character(V1)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


load_dup_set <- function() {
  read.csv(harper_dup_path("Output", "all_paralogs_correct.csv"), stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(FBgn_ID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


gene_universe <- load_gene_universe()
ssdo_ids <- load_ssdo_ids()
sbge_sets <- load_sbge_sets()
retained_ai_ids <- load_mishra_retained_ids()
ai_sets <- load_ai_sets(retained_ai_ids)
dr_set <- load_dr_set()
sSS_set <- load_sss_set()
dup_set <- load_dup_set()


available_feature_sets <- list(
  HFB = unique(sbge_sets$HFB),
  FB = unique(sbge_sets$FB),
  MB = unique(sbge_sets$MB),
  HMB = unique(sbge_sets$HMB),
  ASE = unique(ai_sets$ASE),
  SS_ASE = unique(ai_sets$SS_ASE),
  TS_ASE = unique(ai_sets$TS_ASE),
  DR = unique(dr_set),
  SSS = unique(sSS_set),
  Dup = unique(dup_set)
)


feature_membership <- gene_universe |>
  mutate(IsSSDO = FlyBaseID %in% ssdo_ids)

for (feature_name in names(available_feature_sets)) {
  feature_membership[[feature_name]] <- feature_membership$FlyBaseID %in% available_feature_sets[[feature_name]]
}


feature_sets_long <- bind_rows(lapply(names(available_feature_sets), function(feature_name) {
  data.frame(
    Feature = feature_name,
    FlyBaseID = sort(unique(available_feature_sets[[feature_name]])),
    stringsAsFactors = FALSE
  )
}))


feature_notes <- data.frame(
  Feature = c("HFB", "FB", "MB", "HMB", "ASE", "SS_ASE", "TS_ASE", "DR", "SSS", "Dup"),
  Description = c(
    "Highly female-biased whole-body SBGE genes from Mishra",
    "Female-biased whole-body SBGE genes from Mishra (Moderately female-biased in source file)",
    "Male-biased whole-body SBGE genes from Mishra (Moderately male-biased in source file)",
    "Highly male-biased whole-body SBGE genes from Mishra",
    "Genes with allelic imbalance in the filtered Mishra whole-body sex-reversed AI table (males or females significant)",
    "Genes with significant whole-body ASE in exactly one sex only in the filtered Mishra sex-reversed AI table",
    "Genes with significant ASE in exactly one of heads or gonads in the retained Mishra ASE tables",
    "Genes with dominance-reversed allelic imbalance from Genes.with.Sex.Reversed.AI.tsv",
    "Consistently sex-specific splicing gene set from MaleLimitedEvo / Michelle Liu",
    "Genes with one or more paralogues in the Harper and Morrow duplication analysis"
  ),
  Source = c(
    "Data/inputs/Sex.Biased.Genes.Whole.csv",
    "Data/inputs/Sex.Biased.Genes.Whole.csv",
    "Data/inputs/Sex.Biased.Genes.Whole.csv",
    "Data/inputs/Sex.Biased.Genes.Whole.csv",
    "../Mishra_et_al_Drosophila_ASE/Data_for_R_scripts/Sex.Reversed.AI/Sex.Reversed.AI.test.tsv filtered to ../Genes.for.analysis.final.tsv",
    "../Mishra_et_al_Drosophila_ASE/Data_for_R_scripts/Sex.Reversed.AI/Sex.Reversed.AI.test.tsv filtered to ../Genes.for.analysis.final.tsv",
    "../Mishra_et_al_Drosophila_ASE/Data_for_R_scripts/GLM.outputs/FDR.Sex.Effects.in.AI/{Heads,Gonads}.FDR.test.for.AI.tsv filtered to ../Genes.for.analysis.final.tsv",
    "Data/inputs/Genes.with.Sex.Reversed.AI.tsv",
    "../MaleLimitedEvo/Analyses/DS_analysis/dimorphic.subset.list.txt",
    "Grieshop Lab shared/Data/Harper and Morrow 2025/Output/all_paralogs_correct.csv"
  ),
  Status = c(rep("ready", 10)),
  GeneCount = c(
    length(unique(available_feature_sets$HFB)),
    length(unique(available_feature_sets$FB)),
    length(unique(available_feature_sets$MB)),
    length(unique(available_feature_sets$HMB)),
    length(unique(available_feature_sets$ASE)),
    length(unique(available_feature_sets$SS_ASE)),
    length(unique(available_feature_sets$TS_ASE)),
    length(unique(available_feature_sets$DR)),
    length(unique(available_feature_sets$SSS)),
    length(unique(available_feature_sets$Dup))
  ),
  Notes = c(
    rep(NA_character_, 8),
    "Using dimorphic.subset.list.txt as the final consistently sex-specific splicing (SSS) set from MaleLimitedEvo. The accompanying description indicates that these 2036 genes passed cross-dataset consistency filtering for sex-specific splicing profiles.",
    "Dup is defined here as binary paralogue presence: genes that appear as FBgn_ID in the shared Harper and Morrow all_paralogs_correct.csv table. The current shared table already satisfies the Harper and Morrow paralogue criterion and contains only rows with DIOPT_score >= 2."
  ),
  stringsAsFactors = FALSE
)


feature_sizes <- feature_notes |>
  select(Feature, Status, GeneCount, Description)


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
  feature_sizes,
  summary_path("sa_resolution_feature_sizes.csv"),
  row.names = FALSE
)

write.csv(
  feature_notes,
  summary_path("sa_resolution_feature_notes.csv"),
  row.names = FALSE
)
