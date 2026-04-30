# Compile the core reference tables that define the manuscript's gene classes:
# SSDOs, disease orthologues, all orthologues, and the FlyBase gene reference
# parsed from the current GTF. Most downstream figure scripts depend on these
# outputs.

find_project_root_local <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Plotting")) &&
      dir.exists(file.path(current, "Figures")) &&
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


project_input <- function(..., root = find_project_root_local()) {
  file.path(root, "Data", "inputs", ...)
}


project_output <- function(..., root = find_project_root_local()) {
  file.path(root, "Data", "outputs", ...)
}


project_summary <- function(..., root = find_project_root_local()) {
  file.path(root, "Data", "Summary_tables", ...)
}


# FlyBase-style GTF attributes are stored as semicolon-delimited key/value
# strings; this helper pulls out the specific field we need.
extract_gtf_attribute <- function(attribute_text, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  matches <- regexec(pattern, attribute_text)
  values <- regmatches(attribute_text, matches)
  vapply(values, function(match) {
    if (length(match) >= 2) {
      match[2]
    } else {
      NA_character_
    }
  }, character(1))
}


# Build clean gene-level and transcript-level reference tables from the
# Ensembl/FlyBase GTF. The gene table underlies almost every downstream join.
read_gtf_reference <- function(gtf_path) {
  gtf <- read.delim(
    gtf_path,
    comment.char = "#",
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )

  colnames(gtf) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )

  gene_reference <- gtf |>
    dplyr::filter(feature == "gene") |>
    dplyr::transmute(
      FlyBaseID = extract_gtf_attribute(attribute, "gene_id"),
      Fly.Symbol = extract_gtf_attribute(attribute, "gene_name"),
      gene_biotype = extract_gtf_attribute(attribute, "gene_biotype"),
      chromosome = seqname,
      start = start,
      end = end,
      strand = strand
    ) |>
    dplyr::filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
    dplyr::distinct(FlyBaseID, .keep_all = TRUE) |>
    dplyr::arrange(FlyBaseID)

  transcript_reference <- gtf |>
    dplyr::filter(feature == "transcript") |>
    dplyr::transmute(
      FlyBaseID = extract_gtf_attribute(attribute, "gene_id"),
      transcript_id = extract_gtf_attribute(attribute, "transcript_id"),
      Fly.Symbol = extract_gtf_attribute(attribute, "gene_name"),
      transcript_name = extract_gtf_attribute(attribute, "transcript_name"),
      gene_biotype = extract_gtf_attribute(attribute, "gene_biotype"),
      transcript_biotype = extract_gtf_attribute(attribute, "transcript_biotype"),
      chromosome = seqname
    ) |>
    dplyr::filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
    dplyr::distinct(FlyBaseID, transcript_id, .keep_all = TRUE) |>
    dplyr::arrange(FlyBaseID, transcript_id)

  list(
    gene_reference = gene_reference,
    transcript_reference = transcript_reference
  )
}


# SSDOs combine direct FlyBase disease searches saved as TXT files with
# literature-derived human disease genes mapped to Drosophila via DIOPT.
read_ssdo_text_sources <- function(root) {
  text_files <- list.files(
    project_input("ssdo_inputs", "sex_specific_diseases", root = root),
    pattern = "\\.txt$",
    full.names = TRUE
  )

  dplyr::bind_rows(lapply(text_files, function(path) {
    ids <- readLines(path, warn = FALSE)
    ids <- trimws(ids)
    ids <- ids[nzchar(ids)]

    data.frame(
      FlyBaseID = ids,
      disease = tools::file_path_sans_ext(basename(path)),
      Human.Symbol = NA_character_,
      Fly.Symbol = NA_character_,
      source_type = "FlyBase disease search",
      source_file = basename(path),
      stringsAsFactors = FALSE
    )
  }))
}


# Apply the manuscript's high-confidence DIOPT filter before folding
# literature-derived orthologues into the SSDO set.
read_diopt_source <- function(path, disease_name) {
  data <- read.csv(path, stringsAsFactors = FALSE)
  col_names <- names(data)

  data |>
    dplyr::transmute(
      FlyBaseID = as.character(.data[[col_names[8]]]),
      Human.Symbol = as.character(.data[[col_names[5]]]),
      Fly.Symbol = as.character(.data[[col_names[9]]]),
      DIOPT.Score = as.numeric(.data[[col_names[11]]]),
      Rank = as.character(.data[[col_names[13]]]),
      Best.Score = as.character(.data[[col_names[14]]]),
      Best.Score.Reverse = as.character(.data[[col_names[15]]])
    ) |>
    dplyr::filter(
      Rank == "high",
      Best.Score == "Yes",
      Best.Score.Reverse == "Yes",
      DIOPT.Score >= 3
    ) |>
    dplyr::transmute(
      FlyBaseID,
      disease = disease_name,
      Human.Symbol,
      Fly.Symbol,
      source_type = "Literature orthologue",
      source_file = basename(path)
    )
}


# FlyBase-derived CSVs are not perfectly standardised, so detect the most
# plausible FBgn column rather than assuming a fixed header.
read_fbgn_column <- function(path) {
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  fbgn_counts <- vapply(data, function(column) {
    sum(grepl("^FBgn", column))
  }, integer(1))

  if (max(fbgn_counts) == 0) {
    stop("Could not find a FlyBase gene ID column in ", basename(path), call. = FALSE)
  }

  fbgn_col <- names(fbgn_counts)[which.max(fbgn_counts)]

  data |>
    dplyr::transmute(FlyBaseID = as.character(.data[[fbgn_col]])) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::distinct()
}


# Main wrangling workflow: build the SSDO table, then layer in the broader
# disease-orthologue and orthologue comparators on the same protein-coding gene
# universe.
compile_ssdo_reference_data <- function(project_root = find_project_root_local()) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Missing required R package: dplyr", call. = FALSE)
  }

  dir.create(project_output(root = project_root), recursive = TRUE, showWarnings = FALSE)
  dir.create(project_summary(root = project_root), recursive = TRUE, showWarnings = FALSE)

  gtf_reference <- read_gtf_reference(
    project_input("Drosophila_melanogaster.BDGP6.54.115.gtf", root = project_root)
  )

  text_sources <- read_ssdo_text_sources(project_root) |>
    dplyr::left_join(
      gtf_reference$gene_reference |> dplyr::select(FlyBaseID, Fly.Symbol),
      by = "FlyBaseID",
      suffix = c("", ".gtf")
    ) |>
    dplyr::mutate(Fly.Symbol = dplyr::coalesce(Fly.Symbol, Fly.Symbol.gtf)) |>
    dplyr::select(-Fly.Symbol.gtf)

  literature_sources <- dplyr::bind_rows(
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "diopt_results_PCOS_genes.csv", root = project_root), "PCOS"),
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "diopt_results.csv", root = project_root), "Human_SA"),
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "coronary_artery_disease_orthologs.csv", root = project_root), "CAD"),
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "osteoporosis_orthologs.csv", root = project_root), "Osteoporosis"),
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "alzheimer_disease_orthologs.csv", root = project_root), "Alzheimer"),
    read_diopt_source(project_input("ssdo_inputs", "orthology_sources", "parkinsons_orthologs.csv", root = project_root), "Parkinsons")
  )

  protein_coding_gene_ids <- gtf_reference$gene_reference |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding") |>
    dplyr::pull(FlyBaseID) |>
    unique()

  sexspecific_disease_gene_data <- dplyr::bind_rows(text_sources, literature_sources) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::left_join(
      gtf_reference$gene_reference |> dplyr::select(FlyBaseID, Fly.Symbol.gtf = Fly.Symbol),
      by = "FlyBaseID"
    ) |>
    dplyr::mutate(Fly.Symbol = dplyr::coalesce(Fly.Symbol, Fly.Symbol.gtf)) |>
    dplyr::select(FlyBaseID, Fly.Symbol, Human.Symbol, disease, source_type, source_file) |>
    dplyr::filter(FlyBaseID %in% protein_coding_gene_ids) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID, disease, source_type, source_file)

  ssdo_ids <- sexspecific_disease_gene_data |>
    dplyr::transmute(FlyBaseID) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  raw_disease_ids <- read_fbgn_column(project_input("all_disease_orthologs.csv", root = project_root)) |>
    dplyr::filter(FlyBaseID %in% protein_coding_gene_ids)
  raw_ortholog_ids <- read_fbgn_column(project_input("all_orthologs.csv", root = project_root)) |>
    dplyr::filter(FlyBaseID %in% protein_coding_gene_ids)

  disease_ids <- dplyr::bind_rows(raw_disease_ids, ssdo_ids) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  ortholog_ids <- dplyr::bind_rows(raw_ortholog_ids, disease_ids) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  utils::write.csv(
    sexspecific_disease_gene_data,
    project_output("sexspecific_disease_gene_data.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    ssdo_ids,
    project_output("sexspecific_disease_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    disease_ids,
    project_output("all_disease_ortholog_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    ortholog_ids,
    project_output("all_ortholog_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    gtf_reference$gene_reference,
    project_output("dmel_gene_reference.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    gtf_reference$transcript_reference,
    project_output("dmel_transcript_reference.csv", root = project_root),
    row.names = FALSE
  )

  message("Compiled SSDO rows: ", nrow(sexspecific_disease_gene_data))
  message("Unique SSDO FlyBase IDs: ", nrow(ssdo_ids))
  message("Canonical disease orthologue IDs: ", nrow(disease_ids))
  message("Canonical orthologue IDs: ", nrow(ortholog_ids))

  invisible(
    list(
      sexspecific_disease_gene_data = sexspecific_disease_gene_data,
      ssdo_ids = ssdo_ids,
      disease_ids = disease_ids,
      ortholog_ids = ortholog_ids
    )
  )
}


compile_ssdo_reference_data()
