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
        "Could not locate the project root. Run this script from the Final_Year_Project_github repository.",
        call. = FALSE
      )
    }

    current <- parent
  }
}


project_data <- function(..., root = find_project_root_local()) {
  file.path(root, "Data", ...)
}


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


read_ssdo_text_sources <- function(root) {
  text_files <- list.files(
    project_data("ssdo_inputs", "sex_specific_diseases", root = root),
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


compile_ssdo_reference_data <- function(project_root = find_project_root_local()) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Missing required R package: dplyr", call. = FALSE)
  }

  gtf_reference <- read_gtf_reference(
    project_data("Drosophila_melanogaster.BDGP6.54.115.gtf", root = project_root)
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
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "diopt_results_PCOS_genes.csv", root = project_root), "PCOS"),
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "diopt_results.csv", root = project_root), "Human_SA"),
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "coronary_artery_disease_orthologs.csv", root = project_root), "CAD"),
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "osteoporosis_orthologs.csv", root = project_root), "Osteoporosis"),
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "alzheimer_disease_orthologs.csv", root = project_root), "Alzheimer"),
    read_diopt_source(project_data("ssdo_inputs", "orthology_sources", "parkinsons_orthologs.csv", root = project_root), "Parkinsons")
  )

  sexspecific_disease_gene_data <- dplyr::bind_rows(text_sources, literature_sources) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::left_join(
      gtf_reference$gene_reference |> dplyr::select(FlyBaseID, Fly.Symbol.gtf = Fly.Symbol),
      by = "FlyBaseID"
    ) |>
    dplyr::mutate(Fly.Symbol = dplyr::coalesce(Fly.Symbol, Fly.Symbol.gtf)) |>
    dplyr::select(FlyBaseID, Fly.Symbol, Human.Symbol, disease, source_type, source_file) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID, disease, source_type, source_file)

  ssdo_ids <- sexspecific_disease_gene_data |>
    dplyr::transmute(FlyBaseID) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  raw_disease_ids <- read_fbgn_column(project_data("all_disease_orthologs.csv", root = project_root))
  raw_ortholog_ids <- read_fbgn_column(project_data("all_orthologs.csv", root = project_root))

  disease_ids <- dplyr::bind_rows(raw_disease_ids, ssdo_ids) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  ortholog_ids <- dplyr::bind_rows(raw_ortholog_ids, disease_ids) |>
    dplyr::distinct() |>
    dplyr::arrange(FlyBaseID)

  raw_subset_checks <- data.frame(
    check = c(
      "SSDOs subset raw disease orthologues",
      "Raw disease orthologues subset raw orthologues"
    ),
    passed = c(
      all(ssdo_ids$FlyBaseID %in% raw_disease_ids$FlyBaseID),
      all(raw_disease_ids$FlyBaseID %in% raw_ortholog_ids$FlyBaseID)
    ),
    stringsAsFactors = FALSE
  )

  harmonized_subset_checks <- data.frame(
    check = c(
      "SSDOs subset canonical disease orthologues",
      "Canonical disease orthologues subset canonical orthologues"
    ),
    passed = c(
      all(ssdo_ids$FlyBaseID %in% disease_ids$FlyBaseID),
      all(disease_ids$FlyBaseID %in% ortholog_ids$FlyBaseID)
    ),
    stringsAsFactors = FALSE
  )

  raw_ssdo_missing <- ssdo_ids |>
    dplyr::filter(!(FlyBaseID %in% raw_disease_ids$FlyBaseID))

  raw_disease_missing <- raw_disease_ids |>
    dplyr::filter(!(FlyBaseID %in% raw_ortholog_ids$FlyBaseID))

  utils::write.csv(
    sexspecific_disease_gene_data,
    project_data("sexspecific_disease_gene_data.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    ssdo_ids,
    project_data("sexspecific_disease_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    disease_ids,
    project_data("all_disease_ortholog_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    ortholog_ids,
    project_data("all_ortholog_ids.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    gtf_reference$gene_reference,
    project_data("dmel_gene_reference.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    gtf_reference$transcript_reference,
    project_data("dmel_transcript_reference.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    raw_subset_checks,
    project_data("reference_subset_checks_raw.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    harmonized_subset_checks,
    project_data("reference_subset_checks_harmonized.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    raw_ssdo_missing,
    project_data("ssdo_not_in_raw_disease_orthologs.csv", root = project_root),
    row.names = FALSE
  )
  utils::write.csv(
    raw_disease_missing,
    project_data("raw_disease_not_in_raw_orthologs.csv", root = project_root),
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
      ortholog_ids = ortholog_ids,
      raw_subset_checks = raw_subset_checks,
      harmonized_subset_checks = harmonized_subset_checks
    )
  )
}


compile_ssdo_reference_data()
