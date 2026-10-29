# Build Figure S2 from SBGE-stratified Ruzicka et al. SA-effect summaries.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Sup_plotting_and_analyses")) &&
      dir.exists(file.path(current, "Sup_figures")) &&
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


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})


dir.create(file.path(project_root, "Sup_figures"), recursive = TRUE, showWarnings = FALSE)


ordered_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
# Panels are ordered to mirror the biological progression from all genes to
# increasingly female-biased through increasingly male-biased classes.
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


publication_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1),
      axis.line = element_blank(),
      plot.title = element_text(size = base_size + 3, hjust = 0.5),
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = base_size + 2, face = "bold"),
      legend.text = element_text(size = base_size + 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}


# Retain the same bootstrap summary logic used in the source panels so the
# composite supplementary figure stays numerically consistent with them.
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


# Collapse SNP-level evidence to one gene-level value by interval overlap and
# then taking the strongest signal per gene.
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


attach_sbge_panels <- function(metric_data, sbge_whole) {
  panelled <- metric_data |>
    inner_join(sbge_whole, by = "FlyBaseID") |>
    select(FlyBaseID, Value, GeneType, SexBiasPanel)

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


plot_sbge_panel <- function(summary_data, y_label, title = NULL, tag = NULL, show_x_title = FALSE) {
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
      size = 3.1,
      colour = "black",
      vjust = -0.55
    ) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_x_discrete(
      labels = c("All", "Highly FB", "Female-biased", "Unbiased", "Male-biased", "Highly MB"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    geom_vline(xintercept = 1.5, colour = "black", linewidth = 1) +
    geom_vline(xintercept = 2.5:5.5, colour = "grey60", linewidth = 0.45) +
    labs(
      x = if (show_x_title) "Sex Bias Category" else NULL,
      y = y_label,
      title = title,
      tag = tag,
      colour = "Gene Type"
    ) +
    publication_theme() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.title.x = if (show_x_title) element_text(margin = margin(t = 10)) else element_blank()
    )
}


gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")
ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()

sbge_whole <- read.csv(output_path("mishra_whole_body_sbge_bins.csv"), stringsAsFactors = FALSE) |>
  transmute(
    FlyBaseID = as.character(FlyBaseID),
    SexBiasPanel = recode(Sex.Bias.Category, !!!sbge_label_map)
  ) |>
  filter(!is.na(FlyBaseID), nzchar(FlyBaseID), !is.na(SexBiasPanel)) |>
  distinct(FlyBaseID, .keep_all = TRUE) |>
  mutate(SexBiasPanel = factor(SexBiasPanel, levels = ordered_sbge_panels))


sa_summary <- build_ruzicka_magnitude_values(protein_coding_gene_reference) |>
  assign_gene_types(ssdo_ids = ssdo_ids, disease_ids = disease_ids, ortholog_ids = ortholog_ids) |>
  attach_sbge_panels(sbge_whole = sbge_whole) |>
  summarise_sbge_panels(seed = 301)

figure_s2 <- plot_sbge_panel(
  summary_data = sa_summary,
  y_label = "SA effect",
  show_x_title = TRUE
)

ggsave(
  filename = sup_figure_path("FigS2_sa_mag.png"),
  plot = figure_s2,
  width = 12,
  height = 4.8,
  dpi = 300,
  bg = "white"
)
