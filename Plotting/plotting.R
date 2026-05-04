# Metric plotting workflow
#
# This script is deliberately self-contained so that the data loading,
# four-group assignment, summary-table generation, permutation testing,
# and ggplot figure code all live in one place.
#
# Run after:
#   Rscript Data_wrangling/compile_ssdo_reference_data.R

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


required_files <- c(
  output_path("sexspecific_disease_ids.csv"),
  output_path("all_disease_ortholog_ids.csv"),
  output_path("all_ortholog_ids.csv"),
  output_path("dmel_gene_reference.csv"),
  input_path("FileS7.txt"),
  input_path("SBGEandSSSdataForMBE.csv")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required reference files. Run Data_wrangling/compile_ssdo_reference_data.R first.\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


suppressPackageStartupMessages({
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


# Bootstrap the group mean so the dot/whisker plots show an empirical 95% CI
# even when the underlying metric is non-normal.
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
      Mean = mean(values),
      SD = stats::sd(values),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(summaries) |>
    mutate(GeneType = factor(GeneType, levels = ordered_gene_types)) |>
    arrange(GeneType)
}


# Compare two groups by permuting group labels while holding the observed metric
# values fixed.
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


pairwise_tests_within_metric <- function(data, value_col = "Value", group_col = "GeneType") {
  available_groups <- ordered_gene_types[ordered_gene_types %in% unique(as.character(data[[group_col]]))]
  group_pairs <- combn(available_groups, 2, simplify = FALSE)

  tests <- lapply(group_pairs, function(group_pair) {
    group_1 <- group_pair[1]
    group_2 <- group_pair[2]
    group_1_values <- data[data[[group_col]] == group_1, value_col, drop = TRUE]
    group_2_values <- data[data[[group_col]] == group_2, value_col, drop = TRUE]

    test_result <- two_sample_permutation_test(group_1_values, group_2_values)
    test_result$group_1 <- group_1
    test_result$group_2 <- group_2
    test_result
  })

  bind_rows(tests) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) |>
    select(group_1, group_2, everything())
}


publication_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.line = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 12)),
      plot.margin = margin(12, 12, 12, 12)
    )
}


plot_metric_distribution <- function(summary_data, y_label) {
  ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(
      aes(y = q95, label = SampleSize),
      size = 3.8,
      colour = "black",
      vjust = -0.75
    ) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = "Gene Type", y = y_label, title = NULL) +
    publication_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}


# ---- Main plotting workflow -----------------------------------------------
# Load the metric source tables and canonical gene-set definitions, then export
# the plot-data / permutation-test tables that feed both the standalone panels
# and the assembled manuscript figures.

# Load the metric source tables.
fra <- read.delim(input_path("FileS7.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(fra)[2] <- "FlyBaseID"
fra[sapply(fra, is.infinite)] <- NA
fra <- distinct(fra, FlyBaseID, .keep_all = TRUE)

singh <- read.csv(input_path("SBGEandSSSdataForMBE.csv"), stringsAsFactors = FALSE)
colnames(singh)[2] <- "FlyBaseID"
singh <- distinct(singh, FlyBaseID, .keep_all = TRUE)


# Load the canonical gene-set definitions built by Data_wrangling/.
gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_ids <- gene_reference |>
  transmute(FlyBaseID = as.character(FlyBaseID), gene_biotype = as.character(gene_biotype)) |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding") |>
  pull(FlyBaseID) |>
  unique()

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  transmute(FlyBaseID = as.character(FlyBaseID)) |>
  distinct()

all_disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  transmute(FlyBaseID = as.character(FlyBaseID)) |>
  distinct()

all_ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  transmute(FlyBaseID = as.character(FlyBaseID)) |>
  distinct()

ssdo_gene_ids <- unique(ssdo_ids$FlyBaseID)
disease_only_gene_ids <- setdiff(unique(all_disease_ids$FlyBaseID), ssdo_gene_ids)
ortholog_only_gene_ids <- setdiff(
  unique(all_ortholog_ids$FlyBaseID),
  union(ssdo_gene_ids, unique(all_disease_ids$FlyBaseID))
)


build_metric_data <- function(metric_data, value_col) {
  metric_data |>
    transmute(
      FlyBaseID = as.character(FlyBaseID),
      Value = as.numeric(.data[[value_col]])
    ) |>
    filter(!is.na(FlyBaseID), nzchar(FlyBaseID), FlyBaseID %in% protein_coding_gene_ids) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    filter(!is.na(Value), is.finite(Value)) |>
    mutate(
      GeneType = case_when(
        FlyBaseID %in% ssdo_gene_ids ~ "SSDOs",
        FlyBaseID %in% disease_only_gene_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_gene_ids ~ "Orthologues",
        TRUE ~ "Background"
      ),
      GeneType = factor(GeneType, levels = ordered_gene_types)
    ) |>
    arrange(GeneType, FlyBaseID)
}


metric_specs <- list(
  list(name = "tajD", data = singh, value_col = "tajD.S", y_label = expression("Tajima's " * italic(D))),
  list(name = "rMF", data = singh, value_col = "rmf", y_label = expression(italic(r)[MF])),
  list(name = "rlnd", data = singh, value_col = "piNfract", y_label = expression(tilde(pi)[italic(r) * N])),
  list(name = "DoS", data = fra, value_col = "DoS", y_label = expression(italic(DoS))),
  list(name = "Tissue Specificity", data = singh, value_col = "tao.FML.FA2", y_label = expression("Tissue Specificity (" * tau * ")"))
)

for (metric_spec in metric_specs) {
  raw_metric_data <- build_metric_data(metric_spec$data, metric_spec$value_col)
  summary_metric_data <- bootstrap_group_summaries(raw_metric_data)
  tests_metric_data <- pairwise_tests_within_metric(raw_metric_data)
  metric_plot <- plot_metric_distribution(summary_metric_data, metric_spec$y_label)

  write.csv(
    summary_metric_data,
    summary_path(paste0(metric_spec$name, " Plot Data.csv")),
    row.names = FALSE
  )
  write.csv(
    tests_metric_data,
    summary_path(paste0(metric_spec$name, " Permutation Tests.csv")),
    row.names = FALSE
  )

}
