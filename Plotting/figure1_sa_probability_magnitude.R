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
figure_path <- function(...) file.path(project_root, "Figures", ...)


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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


build_ruzicka_probability_values <- function(gene_reference) {
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


build_neal_probability_values <- function(gene_reference) {
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


build_ruzicka_magnitude_values <- function(gene_reference) {
  assoc <- read.delim(
    input_path("Assoc_Data.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair),
      Value = as.numeric(Wald_Stat)
    ) |>
    filter(!is.na(Value), is.finite(Value))

  map_snp_metric_to_genes(
    snp_values = assoc[, c("chromosome", "position", "Value")],
    gene_reference = gene_reference,
    source_label = "Ruzicka et al. from Assoc_Data.txt (max Wald_Stat)"
  )
}


build_neal_magnitude_values <- function(gene_reference) {
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
      Value = abs(as.numeric(LogOdds))
    ) |>
    filter(!is.na(Value), is.finite(Value))

  map_snp_metric_to_genes(
    snp_values = snps,
    gene_reference = gene_reference,
    source_label = "Neal_Truter_Grieshop from true_cmh_results.txt (max abs LogOdds)"
  )
}


plot_metric_panel <- function(
  summary_data,
  y_label = NULL,
  show_x_axis = TRUE,
  show_y_title = TRUE,
  title = NULL,
  tag = NULL
) {
  y_max <- max(summary_data$q95, na.rm = TRUE)

  p <- ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(aes(y = q95, label = SampleSize), size = 3.9, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(
      x = NULL,
      y = if (show_y_title) y_label else NULL,
      title = title,
      tag = tag
    ) +
    publication_theme()

  if (show_x_axis) {
    p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  } else {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  }

  if (!show_y_title) {
    p <- p + theme(axis.title.y = element_blank())
  }

  p
}


gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")
ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()

ruz_prob_summary <- build_ruzicka_probability_values(protein_coding_gene_reference) |>
  assign_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  bootstrap_group_summaries()

neal_prob_summary <- build_neal_probability_values(protein_coding_gene_reference) |>
  assign_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  bootstrap_group_summaries()

ruz_mag_summary <- build_ruzicka_magnitude_values(protein_coding_gene_reference) |>
  assign_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  bootstrap_group_summaries()

neal_mag_summary <- build_neal_magnitude_values(protein_coding_gene_reference) |>
  assign_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  bootstrap_group_summaries()

panel_a <- plot_metric_panel(
  summary_data = ruz_prob_summary,
  y_label = "Prob. of SA polymorphism",
  show_x_axis = FALSE,
  show_y_title = TRUE,
  title = expression("Ruzicka " * italic("et al.") * " (2019)"),
  tag = "A"
)

panel_b <- plot_metric_panel(
  summary_data = neal_prob_summary,
  y_label = "Prob. of SA polymorphism",
  show_x_axis = FALSE,
  show_y_title = FALSE,
  title = expression("Neal, Truter, " * italic("et al.") * " (unpublished)"),
  tag = "B"
)

panel_c <- plot_metric_panel(
  summary_data = ruz_mag_summary,
  y_label = "SA effect",
  show_x_axis = TRUE,
  show_y_title = TRUE,
  title = NULL,
  tag = "C"
)

panel_d <- plot_metric_panel(
  summary_data = neal_mag_summary,
  y_label = "SA effect",
  show_x_axis = TRUE,
  show_y_title = FALSE,
  title = NULL,
  tag = "D"
)

figure_1 <- (panel_a | panel_b) / (panel_c | panel_d)

ggsave(
  filename = figure_path("Fig1_SA_Probability_and_Magnitude.png"),
  plot = figure_1,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)
