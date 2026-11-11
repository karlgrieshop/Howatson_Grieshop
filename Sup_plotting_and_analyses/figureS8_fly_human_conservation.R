#!/usr/bin/env Rscript

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "fly_human_conservation", "README.md")) &&
        dir.exists(file.path(current, "Data", "outputs"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Run this script from inside the repository.", call. = FALSE)
    }
    current <- parent
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_root <- find_project_root()
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || arguments[[1]] != "--alignment-dir") {
  stop(
    paste(
      "Usage: Rscript Sup_plotting_and_analyses/figureS8_fly_human_conservation.R",
      "--alignment-dir fly_human_conservation/alignment_<job-id>"
    ),
    call. = FALSE
  )
}
alignment_dir <- arguments[[2]]
if (!grepl("^(/|[A-Za-z]:[/\\\\])", alignment_dir)) {
  alignment_dir <- file.path(project_root, alignment_dir)
}
alignment_dir <- normalizePath(alignment_dir, winslash = "/", mustWork = TRUE)
alignment_table_dir <- file.path(alignment_dir, "tables")
summary_dir <- file.path(project_root, "Data", "Summary_tables")
figure_dir <- file.path(project_root, "Sup_figures")
side_figure_dir <- file.path(project_root, "side", "fly_human_conservation")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(side_figure_dir, showWarnings = FALSE, recursive = TRUE)

category_order <- c(
  "SSDOs",
  "Disease orthologues excluding SSDOs",
  "Orthologues excluding disease genes"
)
category_labels <- c(
  "SSDOs",
  "Disease orthologues\n(excluding SSDOs)",
  "Orthologues\n(excluding disease genes)"
)
category_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues excluding SSDOs" = "#0072B2",
  "Orthologues excluding disease genes" = "#009E73"
)
metric_order <- c("identity_percent", "fly_coverage_percent", "human_coverage_percent")
metric_labels <- c(
  identity_percent = "Whole-alignment identity",
  fly_coverage_percent = "Fly protein coverage",
  human_coverage_percent = "Human protein coverage"
)
group_pairs <- list(
  c(category_order[1], category_order[2]),
  c(category_order[1], category_order[3]),
  c(category_order[2], category_order[3])
)
n_boot <- 5000L
n_perm <- 10000L
base_seed <- 1L

main_values <- read.csv(
  file.path(alignment_table_dir, "gene_conservation.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
one_to_one_values <- read.csv(
  file.path(alignment_table_dir, "one_to_one_conservation.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
one_to_one_flag <- tolower(trimws(as.character(one_to_one_values$one_to_one)))
if (any(!one_to_one_flag %in% c("true", "false"))) {
  stop("One-to-one flags must be TRUE/FALSE.", call. = FALSE)
}
one_to_one_values$one_to_one <- one_to_one_flag == "true"
category_counts <- read.csv(
  file.path(alignment_table_dir, "category_counts.csv"),
  stringsAsFactors = FALSE
)

validate_values <- function(data, label) {
  required <- c(
    "FlyBaseID", "category", "human_gene_id", "human_protein_id",
    "one_to_one", metric_order
  )
  if (!all(required %in% names(data))) {
    stop(label, " is missing required columns.", call. = FALSE)
  }
  if (anyDuplicated(data$FlyBaseID) ||
      !setequal(unique(data$category), category_order) ||
      any(!is.finite(as.matrix(data[metric_order]))) ||
      any(as.matrix(data[metric_order]) < 0 | as.matrix(data[metric_order]) > 100)) {
    stop(label, " has duplicated IDs, unexpected categories, or invalid percentages.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_values(main_values, "Main gene table")
validate_values(one_to_one_values, "One-to-one gene table")
if (!all(one_to_one_values$one_to_one %in% TRUE) ||
    !all(one_to_one_values$FlyBaseID %in% main_values$FlyBaseID)) {
  stop("One-to-one table is not a subset of the main table.", call. = FALSE)
}
main_match <- match(one_to_one_values$FlyBaseID, main_values$FlyBaseID)
for (metric in metric_order) {
  if (!isTRUE(all.equal(
    one_to_one_values[[metric]],
    main_values[[metric]][main_match],
    check.attributes = FALSE
  ))) {
    stop("One-to-one values differ from the main table for ", metric, ".", call. = FALSE)
  }
}

required_count_columns <- c(
  "category", "eligible_genes", "excluded_before_alignment", "aligned_pairs",
  "identity_measured", "fly_coverage_measured", "human_coverage_measured",
  "one_to_one_measured"
)
one_to_one_count_columns <- c(
  "one_to_one_eligible", "one_to_one_excluded_before_alignment"
)
if (!all(required_count_columns %in% names(category_counts)) ||
    anyDuplicated(category_counts$category) ||
    !setequal(category_counts$category, category_order)) {
  stop("Category-count table is missing required categories or columns.", call. = FALSE)
}
if (xor(
  one_to_one_count_columns[1] %in% names(category_counts),
  one_to_one_count_columns[2] %in% names(category_counts)
)) {
  stop("One-to-one eligibility counts must be provided together.", call. = FALSE)
}
has_one_to_one_counts <- all(one_to_one_count_columns %in% names(category_counts))
category_counts <- category_counts[
  match(category_order, category_counts$category),
  ,
  drop = FALSE
]
observed_main <- as.integer(table(factor(main_values$category, levels = category_order)))
observed_one_to_one <- as.integer(table(factor(
  one_to_one_values$category,
  levels = category_order
)))
if (!identical(observed_main, as.integer(category_counts$aligned_pairs)) ||
    !identical(observed_one_to_one, as.integer(category_counts$one_to_one_measured)) ||
    any(category_counts$eligible_genes !=
      category_counts$excluded_before_alignment + category_counts$aligned_pairs) ||
    (has_one_to_one_counts && any(category_counts$one_to_one_eligible !=
      category_counts$one_to_one_excluded_before_alignment +
        category_counts$one_to_one_measured)) ||
    any(category_counts$identity_measured != category_counts$aligned_pairs) ||
    any(category_counts$fly_coverage_measured != category_counts$aligned_pairs) ||
    any(category_counts$human_coverage_measured != category_counts$aligned_pairs)) {
  stop("Observed gene universes disagree with the category-count table.", call. = FALSE)
}

summarise_groups <- function(data, analysis_set) {
  bind_rows(lapply(seq_along(metric_order), function(metric_index) {
    metric <- metric_order[metric_index]
    bind_rows(lapply(seq_along(category_order), function(group_index) {
      group <- category_order[group_index]
      values <- data[[metric]][data$category == group]
      seed <- base_seed + metric_index * 100L + group_index
      counts <- category_counts[category_counts$category == group, , drop = FALSE]
      is_main_set <- identical(analysis_set, "All mapped orthologues")
      has_eligible_count <- is_main_set ||
        (identical(analysis_set, "One-to-one sensitivity") && has_one_to_one_counts)
      eligible <- if (is_main_set) {
        counts$eligible_genes
      } else if (has_eligible_count) {
        counts$one_to_one_eligible
      } else {
        NA_integer_
      }
      excluded <- if (is_main_set) {
        counts$excluded_before_alignment
      } else if (has_eligible_count) {
        counts$one_to_one_excluded_before_alignment
      } else {
        NA_integer_
      }
      set.seed(seed)
      bootstrap_means <- replicate(
        n_boot,
        mean(sample(values, size = length(values), replace = TRUE))
      )
      data.frame(
        AnalysisSet = analysis_set,
        Metric = metric,
        Category = group,
        N = length(values),
        Eligible = eligible,
        ExcludedBeforeAlignment = excluded,
        MissingPercent = if (has_eligible_count) 100 * excluded / eligible else NA_real_,
        MeanPercent = mean(values),
        CI95Low = unname(quantile(bootstrap_means, 0.025)),
        CI95High = unname(quantile(bootstrap_means, 0.975)),
        BootstrapReplicates = n_boot,
        Seed = seed,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

test_contrast <- function(data, metric, pair, analysis_set, contrast_order, seed) {
  group1 <- data[data$category == pair[1], , drop = FALSE]
  group2 <- data[data$category == pair[2], , drop = FALSE]
  x <- group1[[metric]]
  y <- group2[[metric]]
  observed_difference <- mean(x) - mean(y)
  pooled <- c(x, y)
  set.seed(seed)
  permuted_differences <- replicate(n_perm, {
    shuffled <- sample(pooled, size = length(pooled), replace = FALSE)
    mean(shuffled[seq_along(x)]) - mean(shuffled[-seq_along(x)])
  })
  p_value <- (sum(abs(permuted_differences) >= abs(observed_difference)) + 1) /
    (n_perm + 1)

  data.frame(
    ContrastOrder = contrast_order,
    AnalysisSet = analysis_set,
    Metric = metric,
    Group1 = pair[1],
    Group2 = pair[2],
    NGroup1 = length(x),
    NGroup2 = length(y),
    MeanGroup1 = mean(x),
    MeanGroup2 = mean(y),
    MeanDifference_pp = observed_difference,
    SharedHumanGenes = length(intersect(
      unique(group1$human_gene_id), unique(group2$human_gene_id)
    )),
    SharedHumanProteins = length(intersect(
      unique(group1$human_protein_id), unique(group2$human_protein_id)
    )),
    PValue = p_value,
    Permutations = n_perm,
    Seed = seed,
    Alternative = "two-sided",
    Method = "Pooled gene-label permutation; fixed group sizes; difference in arithmetic means",
    MultipleTestingFamily = paste(
      analysis_set, metric,
      "three pairwise category comparisons; BH adjusted within this family"
    ),
    stringsAsFactors = FALSE
  )
}

run_analysis_set <- function(data, analysis_set) {
  summaries <- summarise_groups(data, analysis_set)
  tests <- bind_rows(lapply(seq_along(metric_order), function(metric_index) {
    metric <- metric_order[metric_index]
    metric_tests <- bind_rows(lapply(seq_along(group_pairs), function(pair_index) {
      test_contrast(
        data,
        metric,
        group_pairs[[pair_index]],
        analysis_set,
        contrast_order = pair_index,
        seed = base_seed + metric_index * 1000L + pair_index
      )
    }))
    metric_tests$PAdjustedBH <- p.adjust(metric_tests$PValue, method = "BH")
    metric_tests
  }))
  list(summaries = summaries, tests = tests)
}

main_results <- run_analysis_set(main_values, "All mapped orthologues")
one_to_one_results <- run_analysis_set(one_to_one_values, "One-to-one sensitivity")
summaries <- bind_rows(main_results$summaries, one_to_one_results$summaries) |>
  mutate(
    Metric = factor(Metric, levels = metric_order),
    Category = factor(Category, levels = category_order)
  ) |>
  arrange(AnalysisSet, Metric, Category)
tests <- bind_rows(main_results$tests, one_to_one_results$tests) |>
  mutate(Metric = factor(Metric, levels = metric_order)) |>
  arrange(AnalysisSet, Metric, ContrastOrder)

write.csv(
  summaries,
  file.path(summary_dir, "FigS8_group_bootstrap_summaries.csv"),
  row.names = FALSE
)
write.csv(
  tests,
  file.path(summary_dir, "FigS8_pairwise_permutation_tests.csv"),
  row.names = FALSE
)
category_coverage <- category_counts |>
  transmute(
    Category = category,
    Eligible = eligible_genes,
    ExcludedBeforeAlignment = excluded_before_alignment,
    Measured = aligned_pairs,
    MissingPercent = 100 * excluded_before_alignment / eligible_genes,
    OneToOneEligible = if (has_one_to_one_counts) one_to_one_eligible else NA_integer_,
    OneToOneExcludedBeforeAlignment = if (has_one_to_one_counts) {
      one_to_one_excluded_before_alignment
    } else {
      NA_integer_
    },
    OneToOneMissingPercent = if (has_one_to_one_counts) {
      100 * one_to_one_excluded_before_alignment / one_to_one_eligible
    } else {
      NA_real_
    },
    OneToOneMeasured = one_to_one_measured
  )
write.csv(
  category_coverage,
  file.path(summary_dir, "FigS8_category_coverage.csv"),
  row.names = FALSE
)

plot_summaries <- function(summary_data, analysis_set, title) {
  plot_data <- summary_data |>
    filter(AnalysisSet == analysis_set) |>
    mutate(
      Metric = factor(Metric, levels = metric_order, labels = unname(metric_labels)),
      NLabel = ifelse(is.na(Eligible), paste0("n=", N), paste0("n=", N, "/", Eligible))
    )
  ggplot(plot_data, aes(Category, MeanPercent, colour = Category)) +
    geom_errorbar(aes(ymin = CI95Low, ymax = CI95High), width = 0.14, linewidth = 0.75) +
    geom_point(size = 2.8) +
    geom_text(
      aes(y = pmin(CI95High + 3, 99), label = NLabel),
      colour = "black",
      size = 3.2,
      vjust = 0
    ) +
    facet_wrap(~Metric, nrow = 1) +
    scale_x_discrete(labels = category_labels) +
    scale_colour_manual(values = category_palette, drop = FALSE) +
    scale_y_continuous(
      name = "Percent",
      breaks = seq(0, 100, by = 20),
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, title = title, colour = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(colour = "black", angle = 20, hjust = 1, vjust = 1),
      axis.text.y = element_text(colour = "black"),
      axis.title = element_text(size = 13),
      strip.background = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      strip.text = element_text(face = "bold", colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.margin = margin(8, 8, 8, 8)
    )
}

main_plot <- plot_summaries(
  summaries,
  "All mapped orthologues",
  "Fly–human protein conservation"
)
one_to_one_plot <- plot_summaries(
  summaries,
  "One-to-one sensitivity",
  "Fly–human protein conservation: one-to-one sensitivity"
)
ggsave(
  file.path(figure_dir, "FigS8_fly_human_conservation.png"),
  main_plot,
  width = 12,
  height = 5.2,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(side_figure_dir, "one_to_one_sensitivity.png"),
  one_to_one_plot,
  width = 12,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

message("Wrote FigS8 and conservation summaries to Sup_figures/ and Data/Summary_tables/.")
