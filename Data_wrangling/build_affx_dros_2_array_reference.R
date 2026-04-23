find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Data_wrangling")) &&
      dir.exists(file.path(current, "Plotting")) &&
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


project_output <- function(..., root = find_project_root()) {
  file.path(root, "Data", "outputs", ...)
}


suppressPackageStartupMessages({
  library(dplyr)
})


get_affx_dros_2 <- function() {
  if (
    requireNamespace("drosophila2.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)
  ) {
    probes <- AnnotationDbi::keys(drosophila2.db::drosophila2.db, keytype = "PROBEID")

    affx <- AnnotationDbi::select(
      drosophila2.db::drosophila2.db,
      keys = probes,
      columns = c("SYMBOL"),
      keytype = "PROBEID"
    ) |>
      transmute(
        Probe = as.character(PROBEID),
        Symbol = as.character(SYMBOL),
        source_method = "drosophila2.db"
      )
  } else if (requireNamespace("biomaRt", quietly = TRUE)) {
    mart <- biomaRt::useEnsembl("genes", dataset = "dmelanogaster_gene_ensembl")

    affx <- biomaRt::getBM(
      attributes = c("affy_drosophila_2", "external_gene_name"),
      mart = mart
    ) |>
      transmute(
        Probe = as.character(affy_drosophila_2),
        Symbol = as.character(external_gene_name),
        source_method = "biomaRt affy_drosophila_2"
      )
  } else {
    stop(
      "Need either drosophila2.db + AnnotationDbi or biomaRt to construct the Affymetrix Drosophila 2.0 background.",
      call. = FALSE
    )
  }

  affx |>
    filter(!is.na(Probe), nzchar(Probe), !is.na(Symbol), nzchar(Symbol)) |>
    distinct(Probe, Symbol, .keep_all = TRUE) |>
    arrange(Probe, Symbol)
}


build_affx_dros_2_array_reference <- function(project_root = find_project_root()) {
  affx <- get_affx_dros_2()
  dir.create(project_output(root = project_root), recursive = TRUE, showWarnings = FALSE)

  write.csv(
    affx,
    project_output("Affx_Dros_2.0.csv", root = project_root),
    row.names = FALSE
  )

  message("Affymetrix Drosophila 2.0 probe-symbol rows: ", nrow(affx))
  message("Unique symbols: ", dplyr::n_distinct(affx$Symbol))
}


build_affx_dros_2_array_reference()
