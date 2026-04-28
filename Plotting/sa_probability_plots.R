# SA probability plots
#
# This script restores the four-group SA probability plots from the original
# study-specific inputs now present in `Data/`.
#
# - Ruzicka et al.: SNP-level p-values from `Assoc_Data.txt`, mapped to genes
# - Neal_Truter_Grieshop: SNP-level `-log10(p)` values from
#   `true_cmh_results.txt`, mapped to genes

find_project_root <- function(start = getwd()) {
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


project_root <- find_project_root()
input_path <- function(...) file.path(project_root, "Data", "inputs", ...)
output_path <- function(...) file.path(project_root, "Data", "outputs", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
figure_path <- function(...) file.path(project_root, "Figures", ...)


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})


ordered_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")

gene_type_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)


publication_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.line = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.x = element_text(margin = margin(t = 12)),
      plot.margin = margin(32, 12, 12, 24)
    )
}


bootstrap_group_summaries <- function(data, value_col = "Value", group_col = "GeneType", n_boot = 5000, seed = 1L) {
  set.seed(seed)

  split_data <- split(data[[value_col]], data[[group_col]])

  summaries <- lapply(names(split_data), function(group_name) {
    values <- split_data[[group_name]]
    values <- values[!is.na(values)]

    if (length(values) == 0) {
      return(NULL)
    }

    boot_stats <- replicate(
      n_boot,
      mean(sample(values, size = length(values), replace = TRUE))
    )

    data.frame(
      GeneType = group_name,
      q05 = unname(stats::quantile(boot_stats, 0.025)),
      q50 = unname(stats::quantile(boot_stats, 0.5)),
      q95 = unname(stats::quantile(boot_stats, 0.975)),
      SampleSize = length(values),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(summaries) |>
    mutate(GeneType = factor(GeneType, levels = ordered_gene_types)) |>
    arrange(GeneType)
}


two_sample_permutation_test <- function(x, y, n_perm = 10000, seed = 1L) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(
      data.frame(
        obs_diff = NA_real_,
        n_group1 = length(x),
        n_group2 = length(y),
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  set.seed(seed)

  n_x <- length(x)
  pooled <- c(x, y)
  obs_diff <- mean(x) - mean(y)

  perm_diffs <- replicate(n_perm, {
    shuffled <- sample(pooled, size = length(pooled), replace = FALSE)
    mean(shuffled[seq_len(n_x)]) - mean(shuffled[-seq_len(n_x)])
  })

  data.frame(
    obs_diff = obs_diff,
    n_group1 = length(x),
    n_group2 = length(y),
    p_value = (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1),
    stringsAsFactors = FALSE
  )
}


pairwise_tests_against_reference <- function(data, value_col = "Value", group_col = "GeneType", reference_group = "SSDOs") {
  comparison_groups <- setdiff(unique(as.character(data[[group_col]])), reference_group)

  tests <- lapply(comparison_groups, function(group_name) {
    reference_values <- data[data[[group_col]] == reference_group, value_col, drop = TRUE]
    comparison_values <- data[data[[group_col]] == group_name, value_col, drop = TRUE]

    test_result <- two_sample_permutation_test(reference_values, comparison_values)
    test_result$group_1 <- reference_group
    test_result$group_2 <- group_name
    test_result
  })

  bind_rows(tests) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) |>
    select(group_1, group_2, everything())
}


plot_metric_distribution <- function(summary_data, y_label) {
  ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(aes(y = q95, label = SampleSize), size = 3.8, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = "Gene Type", y = y_label, title = NULL) +
    publication_theme() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.title.y = element_text(hjust = 0.42, margin = margin(r = 12))
    )
}


assign_gene_types <- function(raw_values, ssdo_ids, disease_ids, ortholog_ids) {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))

  raw_values |>
    transmute(FlyBaseID = as.character(FlyBaseID), Value = as.numeric(Value), Source = Source) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), !is.na(Value), is.finite(Value)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    mutate(
      GeneType = case_when(
        FlyBaseID %in% ssdo_ids ~ "SSDOs",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      ),
      GeneType = factor(GeneType, levels = ordered_gene_types)
    ) |>
    arrange(GeneType, FlyBaseID)
}


normalise_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  trimws(x)
}


map_snp_metric_to_genes <- function(snp_values, gene_reference, source_label) {
  genes <- gene_reference |>
    transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(start),
      end = as.integer(end),
      FlyBaseID = as.character(FlyBaseID),
      gene_biotype = as.character(gene_biotype)
    ) |>
    filter(
      !is.na(FlyBaseID),
      grepl("^FBgn", FlyBaseID),
      !is.na(chromosome),
      nzchar(chromosome),
      !is.na(start),
      !is.na(end)
    ) |>
    distinct(FlyBaseID, chromosome, start, end, gene_biotype)

  snps <- snp_values |>
    transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(position),
      end = as.integer(position),
      Value = as.numeric(Value)
    ) |>
    filter(
      !is.na(chromosome),
      nzchar(chromosome),
      !is.na(start),
      !is.na(Value),
      is.finite(Value)
    )

  genes_dt <- as.data.table(genes)
  snps_dt <- as.data.table(snps)

  setkey(genes_dt, chromosome, start, end)
  setkey(snps_dt, chromosome, start, end)

  overlaps <- foverlaps(snps_dt, genes_dt, nomatch = 0L)

  as_tibble(overlaps) |>
    filter(!is.na(FlyBaseID)) |>
    group_by(FlyBaseID) |>
    summarise(
      Value = max(Value, na.rm = TRUE),
      Source = source_label,
      .groups = "drop"
    )
}


build_ruzicka_values <- function(gene_reference) {
  assoc <- read.delim(
    input_path("Assoc_Data.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair),
      p_value = as.numeric(Wald_P)
    ) |>
    filter(!is.na(p_value), p_value > 0, p_value <= 1) |>
    mutate(Value = -log10(pmax(p_value, .Machine$double.xmin)))

  map_snp_metric_to_genes(
    snp_values = assoc[, c("chromosome", "position", "Value")],
    gene_reference = gene_reference,
    source_label = "Ruzicka et al. from Assoc_Data.txt (max -log10 Wald_P)"
  )
}


build_neal_values <- function(gene_reference) {
  snps <- read.delim(
    input_path("true_cmh_results.txt"),
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(snps) <- c(
    "Chromosome",
    "Basepair",
    "Reference",
    paste0("Sample_", 1:28),
    "Log10_pvalue",
    "LogOdds"
  )

  snps <- snps |>
    transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair),
      Value = as.numeric(Log10_pvalue)
    ) |>
    filter(!is.na(Value), is.finite(Value))

  map_snp_metric_to_genes(
    snp_values = snps,
    gene_reference = gene_reference,
    source_label = "Neal_Truter_Grieshop from true_cmh_results.txt (max Log10_pvalue)"
  )
}


gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")
ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()

study_specs <- list(
  list(
    name = "Ruzicka et al.",
    values = build_ruzicka_values(protein_coding_gene_reference),
    y_label = expression("Prob. (" * -italic(log)[10] * "(" * italic(P) * ")) of SA polymorphism")
  ),
  list(
    name = "Neal_Truter_Grieshop",
    values = build_neal_values(protein_coding_gene_reference),
    y_label = expression("Prob. (" * -italic(log)[10] * "(" * italic(P) * ")) of SA polymorphism")
  )
)

for (study_spec in study_specs) {
  raw_data <- assign_gene_types(study_spec$values, ssdo_ids, disease_ids, ortholog_ids)
  summary_data <- bootstrap_group_summaries(raw_data)
  tests_data <- pairwise_tests_against_reference(raw_data)
  plot_object <- plot_metric_distribution(summary_data, study_spec$y_label)

  safe_study_name <- gsub("[^A-Za-z0-9]+", "_", study_spec$name)
  safe_study_name <- gsub("^_|_$", "", safe_study_name)

  write.csv(
    summary_data,
    summary_path(paste0(safe_study_name, "_SA_probability_plot_data.csv")),
    row.names = FALSE
  )
  write.csv(
    tests_data,
    summary_path(paste0(safe_study_name, "_SA_probability_permutation_tests.csv")),
    row.names = FALSE
  )

  invisible(plot_object)
}
