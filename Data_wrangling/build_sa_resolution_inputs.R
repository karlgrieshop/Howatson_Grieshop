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
  input_path("Sex.Reversed.AI.test.tsv"),
  input_path("Genes.with.Sex.Reversed.AI.tsv"),
  male_limited_path("Results", "A.f.geno.DS.results.tsv"),
  male_limited_path("Results", "A.m.geno.DS.results.tsv"),
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

  list(SBGE = sbge$FlyBaseID[sbge$Category != "Unbiased"])
}


load_ai_sets <- function() {
  ai <- read.delim(input_path("Sex.Reversed.AI.test.tsv"), stringsAsFactors = FALSE) |>
    transmute(
      FlyBaseID = as.character(geneID),
      p.val.males = as.numeric(p.val.males),
      p.val.females = as.numeric(p.val.females)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE)

  male_sig <- !is.na(ai$p.val.males) & ai$p.val.males < 0.05
  female_sig <- !is.na(ai$p.val.females) & ai$p.val.females < 0.05

  list(
    Cis = ai$FlyBaseID[male_sig | female_sig],
    SS_Cis = ai$FlyBaseID[xor(male_sig, female_sig)]
  )
}


load_dr_set <- function() {
  read.delim(input_path("Genes.with.Sex.Reversed.AI.tsv"), stringsAsFactors = FALSE) |>
    transmute(FlyBaseID = as.character(geneID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    pull(FlyBaseID)
}


load_ds_hits <- function(path) {
  ds <- read.delim(path, stringsAsFactors = FALSE)
  colnames(ds)[2] <- "FlyBaseID"

  ds |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      geneWisePadj = as.numeric(geneWisePadj)
    ) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), !is.na(geneWisePadj), geneWisePadj < 0.1) |>
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
ai_sets <- load_ai_sets()
dr_set <- load_dr_set()
sSS_set <- load_sss_set()
ds_set <- unique(c(
  load_ds_hits(male_limited_path("Results", "A.f.geno.DS.results.tsv")),
  load_ds_hits(male_limited_path("Results", "A.m.geno.DS.results.tsv"))
))
dup_set <- load_dup_set()


available_feature_sets <- list(
  SBGE = unique(sbge_sets$SBGE),
  Cis = unique(ai_sets$Cis),
  SS_Cis = unique(ai_sets$SS_Cis),
  DR = unique(dr_set),
  SSS = unique(sSS_set),
  DS = unique(ds_set),
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
  Feature = c("SBGE", "Cis", "SS_Cis", "DR", "SSS", "DS", "Dup"),
  Description = c(
    "Whole-body SBGE genes from Mishra excluding unbiased genes",
    "Sex-reversed AI test genes with p.val.males < 0.05 or p.val.females < 0.05",
    "Sex-reversed AI test genes significant in exactly one sex only",
    "Genes with dominance-reversed allelic imbalance from Genes.with.Sex.Reversed.AI.tsv",
    "Requested broader sex-specific splicing gene set from MaleLimitedEvo / Michelle Liu",
    "Genes with geneWisePadj < 0.1 in A.f.geno or A.m.geno DS results",
    "Genes annotated as having paralogues in the Harper and Morrow duplication analysis"
  ),
  Source = c(
    "Data/inputs/Sex.Biased.Genes.Whole.csv",
    "Data/inputs/Sex.Reversed.AI.test.tsv",
    "Data/inputs/Sex.Reversed.AI.test.tsv",
    "Data/inputs/Genes.with.Sex.Reversed.AI.tsv",
    "../MaleLimitedEvo/Analyses/DS_analysis/dimorphic.subset.list.txt",
    "../MaleLimitedEvo/Results/A.f.geno.DS.results.tsv + ../MaleLimitedEvo/Results/A.m.geno.DS.results.tsv",
    "Grieshop Lab shared/Data/Harper and Morrow 2025/Output/all_paralogs_correct.csv"
  ),
  Status = c(rep("ready", 7)),
  GeneCount = c(
    length(unique(available_feature_sets$SBGE)),
    length(unique(available_feature_sets$Cis)),
    length(unique(available_feature_sets$SS_Cis)),
    length(unique(available_feature_sets$DR)),
    length(unique(available_feature_sets$SSS)),
    length(unique(available_feature_sets$DS)),
    length(unique(available_feature_sets$Dup))
  ),
  Notes = c(
    rep(NA_character_, 4),
    "Using dimorphic.subset.list.txt as the best available on-disk SSS set because its README describes it as genes with consistent sexual dimorphism in splicing profile across datasets. This yields 2036 genes, which does not match the expected ~3704 and should be checked later against Michelle Liu's original broader list if needed.",
    "DS is the union of unique FlyBase IDs with geneWisePadj < 0.1 from A.f.geno and A.m.geno. C.m.geno is intentionally excluded.",
    "Dup is defined as genes that appear as FBgn_ID in the shared Harper and Morrow all_paralogs_correct.csv table."
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
