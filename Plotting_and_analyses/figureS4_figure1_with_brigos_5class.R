# Supplementary FigS4: SA effect with Brigos as a 5th class.
# Brigos-first assignment rule:
#   Brigos -> SSDOs -> Disease orthologues -> Orthologues -> Background

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  
  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Plotting_and_analyses")) &&
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
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 5-class settings
# ------------------------------------------------------------
ordered_gene_types <- c(
  "Brigos-Baril et al. (2026)",
  "SSDOs",
  "Disease orthologues",
  "Orthologues",
  "Background"
)

gene_type_palette <- c(
  "Brigos-Baril et al. (2026)" = "#CC79A7",
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)

publication_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 2),
      axis.line = element_blank(),
      plot.title = element_text(size = base_size + 3, hjust = 0.5),
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}

bootstrap_group_summaries <- function(data, value_col = "Value", group_col = "GeneType", n_boot = 5000, seed = 1L) {
  set.seed(seed)
  
  split_data <- split(data[[value_col]], data[[group_col]])
  
  summaries <- lapply(names(split_data), function(group_name) {
    values <- split_data[[group_name]]
    values <- values[!is.na(values)]
    
    if (length(values) == 0) return(NULL)
    
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
    return(data.frame(
      obs_diff = NA_real_, n_group1 = length(x), n_group2 = length(y), p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
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

pairwise_tests_within_study <- function(data, value_col = "Value", group_col = "GeneType") {
  available_groups <- ordered_gene_types[ordered_gene_types %in% unique(as.character(data[[group_col]]))]
  group_pairs <- combn(available_groups, 2, simplify = FALSE)
  
  tests <- lapply(group_pairs, function(group_pair) {
    g1 <- group_pair[1]; g2 <- group_pair[2]
    x <- data[data[[group_col]] == g1, value_col, drop = TRUE]
    y <- data[data[[group_col]] == g2, value_col, drop = TRUE]
    test_result <- two_sample_permutation_test(x, y)
    test_result$group_1 <- g1
    test_result$group_2 <- g2
    test_result
  })
  
  bind_rows(tests) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) |>
    select(group_1, group_2, everything())
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
      !is.na(FlyBaseID), grepl("^FBgn", FlyBaseID),
      !is.na(chromosome), nzchar(chromosome), !is.na(start), !is.na(end)
    ) |>
    distinct(FlyBaseID, chromosome, start, end, gene_biotype)
  
  snps <- snp_values |>
    transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(position),
      end = as.integer(position),
      Value = as.numeric(Value)
    ) |>
    filter(!is.na(chromosome), nzchar(chromosome), !is.na(start), !is.na(Value), is.finite(Value))
  
  genes_dt <- as.data.table(genes)
  snps_dt <- as.data.table(snps)
  setkey(genes_dt, chromosome, start, end)
  setkey(snps_dt, chromosome, start, end)
  
  overlaps <- foverlaps(snps_dt, genes_dt, nomatch = 0L)
  
  as_tibble(overlaps) |>
    filter(!is.na(FlyBaseID)) |>
    group_by(FlyBaseID) |>
    summarise(Value = max(Value, na.rm = TRUE), Source = source_label, .groups = "drop")
}

# Brigos-first assignment rule
assign_gene_types_5class <- function(raw_values, ssdo_ids, disease_ids, ortholog_ids, brigos_ids) {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))
  
  raw_values |>
    transmute(FlyBaseID = as.character(FlyBaseID), Value = as.numeric(Value), Source = Source) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), !is.na(Value), is.finite(Value)) |>
    distinct(FlyBaseID, .keep_all = TRUE) |>
    mutate(
      GeneType = case_when(
        FlyBaseID %in% brigos_ids ~ "Brigos-Baril et al. (2026)",
        FlyBaseID %in% ssdo_ids ~ "SSDOs",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      ),
      GeneType = factor(GeneType, levels = ordered_gene_types)
    ) |>
    arrange(GeneType, FlyBaseID)
}

# Study-specific builders (same logic as your current scripts)
build_ruzicka_magnitude_values <- function(gene_reference) {
  assoc <- read.delim(input_path("Assoc_Data.txt"), stringsAsFactors = FALSE, check.names = FALSE) |>
    transmute(chromosome = as.character(Chromosome), position = as.integer(Basepair), Value = as.numeric(Wald_Stat)) |>
    filter(!is.na(Value), is.finite(Value))
  
  map_snp_metric_to_genes(assoc[, c("chromosome", "position", "Value")], gene_reference,
                          "Ruzicka et al. from Assoc_Data.txt (max Wald_Stat)")
}

plot_metric_panel <- function(summary_data, y_label = NULL, show_x_axis = TRUE, show_y_title = TRUE, title = NULL, tag = NULL) {
  p <- ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(aes(y = q95, label = SampleSize), size = 3.8, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_x_discrete(labels = expression(
      atop("Brigos-Baril", italic("et al.") * " (2026)"),
      "SSDOs", "Disease orthologues", "Orthologues", "Background"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = NULL, y = if (show_y_title) y_label else NULL, title = title, tag = tag) +
    publication_theme()
  
  if (show_x_axis) {
    p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  }
  
  if (!show_y_title) p <- p + theme(axis.title.y = element_blank())
  
  p
}

# ------------------------------------------------------------
# Main run: rebuild summaries + figure
# ------------------------------------------------------------
gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
brigos_ids <- read.csv(output_path("brigos_barril_ids_strict.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()

# Magnitude summaries
ruz_mag_raw <- assign_gene_types_5class(
  build_ruzicka_magnitude_values(protein_coding_gene_reference),
  ssdo_ids, disease_ids, ortholog_ids, brigos_ids
)

ruz_mag_summary <- bootstrap_group_summaries(ruz_mag_raw, seed = 21L)

ruz_mag_tests <- pairwise_tests_within_study(ruz_mag_raw)

write.csv(ruz_mag_summary, summary_path("Ruzicka_et_al_SA_magnitude_plot_data_brigos5.csv"), row.names = FALSE)
write.csv(ruz_mag_tests, summary_path("Ruzicka_et_al_SA_magnitude_permutation_tests_brigos5.csv"), row.names = FALSE)

figure_1_brigos5 <- plot_metric_panel(
  summary_data = ruz_mag_summary,
  y_label = "SA effect",
  show_x_axis = TRUE,
  show_y_title = TRUE,
  title = NULL
)

ggsave(
  filename = sup_figure_path("FigS4_SA_Magnitude_Brigos5.png"),
  plot = figure_1_brigos5,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

message("Done:")
message(" - Rebuilt SA-effect summaries with Brigos 5th class")
message(" - Rebuilt supplementary Figure S4: Sup_figures/FigS4_SA_Magnitude_Brigos5.png")