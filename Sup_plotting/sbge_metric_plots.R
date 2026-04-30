# SBGE-stratified supplementary metric plots. These repeat the main metric
# comparisons within each whole-body SBGE class so that the supplement can show
# whether the four-way contrasts are consistent across sex-bias backgrounds.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Sup_plotting")) &&
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
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", "SBGE_supplement", ...)


required_files <- c(
  output_path("sexspecific_disease_ids.csv"),
  output_path("all_disease_ortholog_ids.csv"),
  output_path("all_ortholog_ids.csv"),
  output_path("dmel_gene_reference.csv"),
  output_path("mishra_whole_body_sbge_bins.csv"),
  input_path("FileS7.txt"),
  input_path("SBGEandSSSdataForMBE.csv")
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required reference files. Run Data_wrangling/compile_ssdo_reference_data.R first and make sure the Mishra SBGE file is present.\n",
      "Run Data_wrangling/build_mishra_whole_body_sbge_bins.R to rebuild the active Mishra whole-body SBGE bins.\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


dir.create(file.path(project_root, "Sup_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "Data", "Summary_tables", "SBGE_supplement"), recursive = TRUE, showWarnings = FALSE)


ordered_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
ordered_sbge_panels <- c(
  "All",
  "Highly female-biased",
  "Female-biased",
  "Unbiased",
  "Male-biased",
  "Highly male-biased"
)

gene_type_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)

sbge_label_map <- c(
  "Highly female-biased" = "Highly female-biased",
  "Moderately female-biased" = "Female-biased",
  "Unbiased" = "Unbiased",
  "Moderately male-biased" = "Male-biased",
  "Highly male-biased" = "Highly male-biased"
)


# Bootstrap the group mean within each SBGE panel to obtain empirical 95% CIs
# for the supplementary point-range summaries.
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


# Permutation tests remain SSDO-versus-control contrasts, but are repeated
# within each SBGE stratum.
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


pairwise_tests_against_reference <- function(data, value_col = "Value", group_col = "GeneType", reference_group = "SSDOs", seed = 1L) {
  comparison_groups <- setdiff(unique(as.character(data[[group_col]])), reference_group)

  tests <- lapply(seq_along(comparison_groups), function(i) {
    group_name <- comparison_groups[[i]]
    reference_values <- data[data[[group_col]] == reference_group, value_col, drop = TRUE]
    comparison_values <- data[data[[group_col]] == group_name, value_col, drop = TRUE]

    test_result <- two_sample_permutation_test(reference_values, comparison_values, seed = seed + i)
    test_result$group_1 <- reference_group
    test_result$group_2 <- group_name
    test_result
  })

  bind_rows(tests) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) |>
    select(group_1, group_2, everything())
}


publication_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.line = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}


plot_metric_distribution_by_sbge <- function(summary_data, y_label) {
  dodge_width <- 0.55

  ggplot(summary_data, aes(x = SexBiasPanel, y = q50, colour = GeneType, group = GeneType)) +
    geom_pointrange(
      aes(ymin = q05, ymax = q95),
      position = position_dodge(width = dodge_width),
      linewidth = 0.9
    ) +
    geom_text(
      aes(y = q95, label = SampleSize),
      position = position_dodge(width = dodge_width),
      size = 2.8,
      colour = "black",
      vjust = -0.55
    ) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_x_discrete(
      labels = c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = "Sex Bias Category", y = y_label, colour = "Gene Type", title = NULL) +
    geom_vline(xintercept = 1.5, colour = "black", linewidth = 1) +
    geom_vline(xintercept = 2.5:5.5, colour = "grey60", linewidth = 0.45) +
    publication_theme() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      strip.background = element_blank()
    )
}


# ---- Main supplementary workflow ------------------------------------------

# Load the metric source tables.
fra <- read.delim(input_path("FileS7.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(fra)[2] <- "FlyBaseID"
fra[sapply(fra, is.infinite)] <- NA
fra <- distinct(fra, FlyBaseID, .keep_all = TRUE)

singh <- read.csv(input_path("SBGEandSSSdataForMBE.csv"), stringsAsFactors = FALSE)
colnames(singh)[2] <- "FlyBaseID"
singh <- distinct(singh, FlyBaseID, .keep_all = TRUE)

sbge_whole <- read.csv(output_path("mishra_whole_body_sbge_bins.csv"), stringsAsFactors = FALSE) |>
  transmute(
    FlyBaseID = as.character(FlyBaseID),
    SexBiasPanel = recode(Sex.Bias.Category, !!!sbge_label_map),
    SBGELog2FoldChange = as.numeric(log2FoldChange)
  ) |>
  filter(!is.na(FlyBaseID), nzchar(FlyBaseID), !is.na(SexBiasPanel)) |>
  distinct(FlyBaseID, .keep_all = TRUE) |>
  mutate(SexBiasPanel = factor(SexBiasPanel, levels = ordered_sbge_panels))


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


attach_sbge_panels <- function(metric_data) {
  panelled <- metric_data |>
    inner_join(sbge_whole, by = "FlyBaseID") |>
    select(FlyBaseID, Value, GeneType, SexBiasPanel, SBGELog2FoldChange)

  bind_rows(
    panelled,
    panelled |>
      mutate(SexBiasPanel = factor("All", levels = ordered_sbge_panels))
  ) |>
    mutate(SexBiasPanel = factor(as.character(SexBiasPanel), levels = ordered_sbge_panels))
}


summarise_sbge_panels <- function(panelled_metric_data, seed = 1L) {
  panel_summaries <- lapply(seq_along(ordered_sbge_panels), function(i) {
    panel_name <- ordered_sbge_panels[[i]]
    panel_data <- panelled_metric_data |>
      filter(SexBiasPanel == panel_name)

    if (nrow(panel_data) == 0) {
      return(NULL)
    }

    bootstrap_group_summaries(panel_data, n_boot = 5000, seed = seed + i) |>
      mutate(SexBiasPanel = factor(panel_name, levels = ordered_sbge_panels))
  })

  bind_rows(panel_summaries) |>
    select(SexBiasPanel, everything())
}


test_sbge_panels <- function(panelled_metric_data, seed = 1L) {
  panel_tests <- lapply(seq_along(ordered_sbge_panels), function(i) {
    panel_name <- ordered_sbge_panels[[i]]
    panel_data <- panelled_metric_data |>
      filter(SexBiasPanel == panel_name)

    if (nrow(panel_data) == 0 || dplyr::n_distinct(panel_data$GeneType) < 2 || !"SSDOs" %in% panel_data$GeneType) {
      return(NULL)
    }

    pairwise_tests_against_reference(panel_data, seed = seed + i) |>
      mutate(SexBiasPanel = factor(panel_name, levels = ordered_sbge_panels))
  })

  bind_rows(panel_tests) |>
    select(SexBiasPanel, everything())
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
  sbge_metric_data <- attach_sbge_panels(raw_metric_data)
  summary_metric_data <- summarise_sbge_panels(sbge_metric_data)
  tests_metric_data <- test_sbge_panels(sbge_metric_data)
  metric_plot <- plot_metric_distribution_by_sbge(summary_metric_data, metric_spec$y_label)

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
  invisible(metric_plot)

}
