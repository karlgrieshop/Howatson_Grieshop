find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Analyses")) &&
      dir.exists(file.path(current, "Figures")) &&
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
data_path <- function(...) file.path(project_root, "Data", ...)
figure_path <- function(...) file.path(project_root, "Figures", ...)


suppressPackageStartupMessages({
  library(dplyr)
})


# This file is a scaffold for future enrichment tests asking whether SSDOs are
# enriched for genes from other studies that show evidence of sexually
# antagonistic resolution. It intentionally does not run any tests yet because
# the specific input files and their universe definitions have not been fixed.

ssdo_ids <- read.csv(data_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()

gene_reference <- read.csv(data_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)


run_resolution_enrichment <- function(universe_ids, candidate_ids, study_name, universe_type, source_note) {
  analysis_df <- data.frame(FlyBaseID = unique(universe_ids), stringsAsFactors = FALSE) |>
    mutate(
      IsSSDO = FlyBaseID %in% ssdo_ids,
      IsCandidate = FlyBaseID %in% candidate_ids
    )

  contingency <- table(analysis_df$IsSSDO, analysis_df$IsCandidate)
  fisher_result <- fisher.test(contingency)

  data.frame(
    Study = study_name,
    UniverseType = universe_type,
    Source = source_note,
    UniverseSize = nrow(analysis_df),
    SSDOCount = sum(analysis_df$IsSSDO),
    CandidateCount = sum(analysis_df$IsCandidate),
    OverlapCount = sum(analysis_df$IsSSDO & analysis_df$IsCandidate),
    OddsRatio = unname(fisher_result$estimate),
    PValue = fisher_result$p.value,
    ConfLow = fisher_result$conf.int[1],
    ConfHigh = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )
}


# Planned structure for future use:
#
# 1. Load a resolution-candidate gene set from another study.
# 2. Define two universes where possible:
#    - broad whole-genome / protein-coding background
#    - study-matched tested-gene background
# 3. Run run_resolution_enrichment() for each study and universe type.
# 4. Write outputs to Data/Summary_tables/ and, if useful, a figure to Figures/.
#
# This script is intentionally non-executing for now so the repository has a
# clean, documented place for the next round of SA-resolution enrichment work.
