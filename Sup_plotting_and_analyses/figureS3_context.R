#!/usr/bin/env Rscript

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Plotting_and_analyses")) &&
      dir.exists(file.path(current, "Sup_plotting_and_analyses")) &&
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
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)
source(file.path(project_root, "Data_wrangling", "load_grieshop_candidate_ids.R"))


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})


dir.create(sup_figure_path(), recursive = TRUE, showWarnings = FALSE)


ordered_gene_types <- c(
  "SSDOs",
  "Innocenti and Morrow SA candidates",
  "Grieshop et al. (2025) SA candidates",
  "Ruzicka et al. SA candidates"
)
original_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")

gene_type_palette <- c(
  "SSDOs" = "#D55E00",
  "Innocenti and Morrow SA candidates" = "#7B3294",
  "Grieshop et al. (2025) SA candidates" = "#333333",
  "Ruzicka et al. SA candidates" = "#E6AB02"
)
reference_line_palette <- c(
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)

gene_type_labels <- expression(
  "SSDOs",
  "Innocenti and Morrow candidates",
  "Grieshop " * italic("et al.") * " candidates",
  "Ruzicka " * italic("et al.") * " candidates"
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


bootstrap_group_summaries <- function(data, value_col = "Value", group_col = "GeneType", group_levels = ordered_gene_types, n_boot = 5000, seed = 1L) {
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
    mutate(GeneType = factor(GeneType, levels = group_levels)) |>
    arrange(GeneType)
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


map_symbols_to_flybase <- function(symbols, gene_reference) {
  symbol_map <- gene_reference |>
    filter(!is.na(Fly.Symbol), nzchar(Fly.Symbol)) |>
    select(FlyBaseID, Fly.Symbol) |>
    distinct() |>
    add_count(Fly.Symbol, name = "n_ids") |>
    filter(n_ids == 1) |>
    select(-n_ids)

  data.frame(Symbol = as.character(symbols), stringsAsFactors = FALSE) |>
    distinct() |>
    left_join(symbol_map, by = c("Symbol" = "Fly.Symbol"))
}


load_innocenti_candidate_ids <- function(gene_reference) {
  read.csv(
    input_path("Innocenti_and_Morrow_ANTAGONISTIC.csv"),
    stringsAsFactors = FALSE
  ) |>
    transmute(Symbol = as.character(Symbol)) |>
    filter(!is.na(Symbol), nzchar(Symbol)) |>
    distinct() |>
    pull(Symbol) |>
    map_symbols_to_flybase(gene_reference = gene_reference) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    pull(FlyBaseID) |>
    unique()
}


load_ruzicka_candidate_ids <- function() {
  read.csv(
    input_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    pull(FlyBaseID) |>
    unique()
}


assign_overlay_gene_types <- function(raw_values, group_ids) {
  bind_rows(lapply(ordered_gene_types, function(group_name) {
    ids <- unique(as.character(group_ids[[group_name]]))

    raw_values |>
      transmute(FlyBaseID = as.character(FlyBaseID), Value = as.numeric(Value), Source = Source) |>
      filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), !is.na(Value), is.finite(Value)) |>
      distinct(FlyBaseID, .keep_all = TRUE) |>
      filter(FlyBaseID %in% ids) |>
      mutate(GeneType = group_name)
  })) |>
    mutate(GeneType = factor(GeneType, levels = ordered_gene_types)) |>
    arrange(GeneType, FlyBaseID)
}


assign_original_gene_types <- function(raw_values, ssdo_ids, disease_ids, ortholog_ids) {
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
      GeneType = factor(GeneType, levels = original_gene_types)
    ) |>
    arrange(GeneType, FlyBaseID)
}


plot_metric_panel <- function(
  summary_data,
  reference_lines = NULL,
  y_label = NULL,
  show_x_axis = TRUE,
  show_y_title = TRUE,
  title = NULL,
  tag = NULL
) {
  y_values <- summary_data$q95
  if (!is.null(reference_lines) && nrow(reference_lines) > 0) {
    y_values <- c(y_values, reference_lines$yintercept)
  }

  ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    {
      if (!is.null(reference_lines) && nrow(reference_lines) > 0) {
        geom_hline(
          data = reference_lines,
          aes(yintercept = yintercept, colour = LineGroup),
          inherit.aes = FALSE,
          linewidth = 0.8
        )
      }
    } +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(aes(y = q95, label = SampleSize), size = 3.9, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = c(gene_type_palette, reference_line_palette), drop = FALSE) +
    scale_x_discrete(labels = gene_type_labels, drop = FALSE) +
    scale_y_continuous(
      expand = expansion(mult = c(0.03, 0.14)),
      limits = c(NA, max(y_values, na.rm = TRUE) * 1.02)
    ) +
    labs(
      x = NULL,
      y = if (show_y_title) y_label else NULL,
      title = title,
      tag = tag
    ) +
    publication_theme() +
    {
      if (show_x_axis) {
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
      } else {
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        )
      }
    } +
    {
      if (!show_y_title) theme(axis.title.y = element_blank()) else theme()
    }
}


gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
innocenti_candidate_ids <- intersect(load_innocenti_candidate_ids(gene_reference), protein_coding_gene_reference$FlyBaseID)
ruzicka_candidate_ids <- intersect(load_ruzicka_candidate_ids(), protein_coding_gene_reference$FlyBaseID)
grieshop_ids <- load_grieshop_candidate_ids(input_path("DE.summary.tsv"))
grieshop_group <- "Grieshop et al. (2025) SA candidates"

group_ids <- list(
  "SSDOs" = ssdo_ids,
  "Innocenti and Morrow SA candidates" = innocenti_candidate_ids,
  "Grieshop et al. (2025) SA candidates" = intersect(grieshop_ids, protein_coding_gene_reference$FlyBaseID),
  "Ruzicka et al. SA candidates" = ruzicka_candidate_ids
)

ruz_mag_values <- build_ruzicka_magnitude_values(protein_coding_gene_reference)

overlay_values <- assign_overlay_gene_types(ruz_mag_values, group_ids = group_ids)
# Keep the original groups' bootstrap draws unchanged when adding a candidate set.
base_groups <- setdiff(ordered_gene_types, grieshop_group)
ruz_mag_summary <- bind_rows(
  overlay_values |>
    filter(GeneType != grieshop_group) |>
    mutate(GeneType = factor(GeneType, levels = base_groups)) |>
    bootstrap_group_summaries(group_levels = base_groups, seed = 13L),
  overlay_values |>
    filter(GeneType == grieshop_group) |>
    droplevels() |>
    bootstrap_group_summaries(group_levels = grieshop_group, seed = 13L)
) |>
  mutate(GeneType = factor(GeneType, levels = ordered_gene_types)) |>
  arrange(GeneType)
ruz_mag_reference <- ruz_mag_values |>
  assign_original_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  bootstrap_group_summaries(group_levels = original_gene_types, seed = 23L) |>
  filter(as.character(GeneType) %in% names(reference_line_palette)) |>
  transmute(LineGroup = as.character(GeneType), yintercept = q50)

grieshop_audit <- data.frame(FlyBaseID = grieshop_ids) |>
  mutate(InProteinCodingUniverse = FlyBaseID %in% protein_coding_gene_reference$FlyBaseID) |>
  left_join(ruz_mag_values |> select(FlyBaseID, Value), by = "FlyBaseID") |>
  mutate(
    IncludedInPlot = InProteinCodingUniverse & !is.na(Value),
    SourceFile = "Data/inputs/DE.summary.tsv",
    Selection = "Sig == TRUE"
  )
dir.create(summary_path(), recursive = TRUE, showWarnings = FALSE)
write.csv(grieshop_audit, summary_path("FigS3_Grieshop_candidate_audit.csv"), row.names = FALSE)
write.csv(overlay_values, summary_path("FigS3_gene_level_values.csv"), row.names = FALSE)
write.csv(ruz_mag_summary, summary_path("FigS3_bootstrap_summaries.csv"), row.names = FALSE)
write.csv(ruz_mag_reference, summary_path("FigS3_reference_lines.csv"), row.names = FALSE)

figure_context <- plot_metric_panel(
  summary_data = ruz_mag_summary,
  reference_lines = ruz_mag_reference,
  y_label = "SA effect",
  show_x_axis = TRUE,
  show_y_title = TRUE,
  title = NULL
)

ggsave(
  filename = sup_figure_path("FigS3_sa_mag_context.png"),
  plot = figure_context,
  width = 11,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
