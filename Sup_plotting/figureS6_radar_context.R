#!/usr/bin/env Rscript

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Analyses")) &&
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
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)

suppressPackageStartupMessages({
  library(dplyr)
})

if (!requireNamespace("fmsb", quietly = TRUE)) {
  stop("Install the 'fmsb' package to draw the supplementary radar chart.", call. = FALSE)
}

dir.create(sup_figure_path(), recursive = TRUE, showWarnings = FALSE)


normalise_membership_columns <- function(data) {
  logical_cols <- setdiff(colnames(data), c("FlyBaseID", "Fly.Symbol", "gene_biotype", "chromosome"))
  for (col_name in logical_cols) {
    data[[col_name]] <- as.logical(data[[col_name]])
  }
  data
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


load_neal_candidate_ids <- function() {
  read.csv(
    input_path("gene_hits_snps_90.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(ensembl_gene_id)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    pull(FlyBaseID) |>
    unique()
}


build_radar_components <- function(data, group_var, group_order, feature_order, feature_labels) {
  radar_summary <- bind_rows(lapply(group_order, function(group_name) {
    subset_df <- data |>
      filter(.data[[group_var]] == group_name)

      data.frame(
        Group = group_name,
        Feature = feature_order,
        FeatureLabel = unname(feature_labels[feature_order]),
        Percent = vapply(feature_order, function(feature_name) {
          100 * mean(subset_df[[feature_name]], na.rm = TRUE)
        }, numeric(1)),
        NConsidered = vapply(feature_order, function(feature_name) {
          sum(!is.na(subset_df[[feature_name]]))
        }, integer(1)),
        stringsAsFactors = FALSE
      )
    }))

  if (any(radar_summary$NConsidered == 0L)) {
    empty_rows <- radar_summary |>
      filter(NConsidered == 0L) |>
      transmute(Label = paste(Group, Feature, sep = " / ")) |>
      pull(Label)

    stop(
      paste(
        "At least one supplementary group-feature combination has zero study-matched genes:",
        paste(empty_rows, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  max_percent <- ceiling(max(radar_summary$Percent, na.rm = TRUE) / 10) * 10
  max_percent <- max(max_percent, 10)
  radial_step <- if (max_percent <= 50) 10 else 20
  radial_breaks <- seq(0, max_percent, by = radial_step)
  plot_group_order <- rev(group_order)

  radar_matrix <- rbind(
    rep(max_percent, length(feature_order)),
    rep(0, length(feature_order)),
    do.call(
      rbind,
      lapply(plot_group_order, function(group_name) {
        radar_summary |>
          filter(Group == group_name) |>
          arrange(match(Feature, feature_order)) |>
          pull(Percent)
      })
    )
  )

  radar_df <- as.data.frame(radar_matrix, stringsAsFactors = FALSE)
  colnames(radar_df) <- unname(feature_labels[feature_order])
  rownames(radar_df) <- c("max", "min", plot_group_order)

  list(
    summary = radar_summary,
    radar_df = radar_df,
    radial_breaks = radial_breaks,
    plot_group_order = plot_group_order
  )
}


build_overlay_radar_components <- function(data, group_ids, group_order, feature_order, feature_labels) {
  radar_summary <- bind_rows(lapply(group_order, function(group_name) {
    subset_ids <- unique(as.character(group_ids[[group_name]]))
    subset_df <- data |>
      filter(FlyBaseID %in% subset_ids)

    if (nrow(subset_df) == 0) {
      stop(sprintf("No genes available for overlay group '%s'.", group_name), call. = FALSE)
    }

      data.frame(
        Group = group_name,
        Feature = feature_order,
        FeatureLabel = unname(feature_labels[feature_order]),
        Percent = vapply(feature_order, function(feature_name) {
          100 * mean(subset_df[[feature_name]], na.rm = TRUE)
        }, numeric(1)),
        NConsidered = vapply(feature_order, function(feature_name) {
          sum(!is.na(subset_df[[feature_name]]))
        }, integer(1)),
        stringsAsFactors = FALSE
      )
    }))

  if (any(radar_summary$NConsidered == 0L)) {
    empty_rows <- radar_summary |>
      filter(NConsidered == 0L) |>
      transmute(Label = paste(Group, Feature, sep = " / ")) |>
      pull(Label)

    stop(
      paste(
        "At least one supplementary overlay group-feature combination has zero study-matched genes:",
        paste(empty_rows, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  max_percent <- ceiling(max(radar_summary$Percent, na.rm = TRUE) / 10) * 10
  max_percent <- max(max_percent, 10)
  radial_step <- if (max_percent <= 50) 10 else 20
  radial_breaks <- seq(0, max_percent, by = radial_step)
  plot_group_order <- rev(group_order)

  radar_matrix <- rbind(
    rep(max_percent, length(feature_order)),
    rep(0, length(feature_order)),
    do.call(
      rbind,
      lapply(plot_group_order, function(group_name) {
        radar_summary |>
          filter(Group == group_name) |>
          arrange(match(Feature, feature_order)) |>
          pull(Percent)
      })
    )
  )

  radar_df <- as.data.frame(radar_matrix, stringsAsFactors = FALSE)
  colnames(radar_df) <- unname(feature_labels[feature_order])
  rownames(radar_df) <- c("max", "min", plot_group_order)

  list(
    summary = radar_summary,
    radar_df = radar_df,
    radial_breaks = radial_breaks,
    plot_group_order = plot_group_order
  )
}


draw_radar_panel <- function(components, panel_title, panel_tag, panel_palette) {
  fmsb::radarchart(
    components$radar_df,
    axistype = 1,
    pcol = panel_palette[components$plot_group_order],
    pfcol = grDevices::adjustcolor(panel_palette[components$plot_group_order], alpha.f = 0.20),
    plwd = 4,
    plty = 1,
    cglcol = "grey78",
    cglty = 1,
    axislabcol = "grey45",
    caxislabels = as.character(components$radial_breaks),
    cglwd = 0.8,
    calcex = 1.35,
    vlcex = 1.25
  )

  if (!is.null(panel_title) && nzchar(panel_title)) {
    title(main = panel_title, line = 0.7, cex.main = 1.3, font.main = 1)
  }

  if (!is.null(panel_tag) && nzchar(panel_tag)) {
    mtext(panel_tag, side = 3, line = 1.4, adj = 0.02, font = 2, cex = 1.25)
  }
}


assign_resolution_groups <- function(data, ssdo_ids, disease_ids, ortholog_ids) {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))

  data |>
    mutate(
      PlotGroup = case_when(
        FlyBaseID %in% ssdo_ids ~ "SSDOs",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      )
    )
}


assign_split_candidate_groups <- function(data, candidate_ids, ssdo_ids, disease_ids, ortholog_ids) {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))
  excluded_ssdo_ids <- setdiff(ssdo_ids, candidate_ids)

  data |>
    filter(!(FlyBaseID %in% excluded_ssdo_ids)) |>
    mutate(
      PlotGroup = case_when(
        FlyBaseID %in% candidate_ids ~ "Candidate SA genes",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      )
    )
}


resolution_membership_path <- output_path("sa_resolution_feature_membership.csv")
resolution_notes_path <- summary_path("sa_resolution_feature_notes.csv")

if (!file.exists(resolution_membership_path) || !file.exists(resolution_notes_path)) {
  stop(
    "Missing SA-resolution input tables. Run Data_wrangling/build_sa_resolution_inputs.R first.",
    call. = FALSE
  )
}

resolution_notes <- read.csv(resolution_notes_path, stringsAsFactors = FALSE)
feature_order <- resolution_notes |>
  filter(Status == "ready") |>
  pull(Feature)

feature_labels <- setNames(feature_order, feature_order)
feature_labels[c("SS_ASE", "TS_ASE")] <- c("SS-ASE", "TS-ASE")

resolution_membership <- read.csv(resolution_membership_path, stringsAsFactors = FALSE) |>
  normalise_membership_columns() |>
  filter(gene_biotype == "protein_coding") |>
  select(FlyBaseID, IsSSDO, all_of(feature_order))

gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()

main_group_order <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
candidate_group_order <- c("Candidate SA genes", "Disease orthologues", "Orthologues", "Background")

main_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)
candidate_palette <- c(
  "Candidate SA genes" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)
overlay_group_order <- c(
  "SSDOs",
  "Innocenti and Morrow SA candidates",
  "Ruzicka et al. SA candidates",
  "Neal, Truter et al. SA candidates"
)
overlay_palette <- c(
  "SSDOs" = "#D55E00",
  "Innocenti and Morrow SA candidates" = "#7B3294",
  "Ruzicka et al. SA candidates" = "#E6AB02",
  "Neal, Truter et al. SA candidates" = "#B2182B"
)

innocenti_candidate_ids <- intersect(load_innocenti_candidate_ids(gene_reference), resolution_membership$FlyBaseID)
ruzicka_candidate_ids <- intersect(load_ruzicka_candidate_ids(), resolution_membership$FlyBaseID)
neal_candidate_ids <- intersect(load_neal_candidate_ids(), resolution_membership$FlyBaseID)

main_plot_data <- assign_resolution_groups(resolution_membership, ssdo_ids, disease_ids, ortholog_ids)
innocenti_plot_data <- assign_split_candidate_groups(resolution_membership, innocenti_candidate_ids, ssdo_ids, disease_ids, ortholog_ids)
ruzicka_plot_data <- assign_split_candidate_groups(resolution_membership, ruzicka_candidate_ids, ssdo_ids, disease_ids, ortholog_ids)
neal_plot_data <- assign_split_candidate_groups(resolution_membership, neal_candidate_ids, ssdo_ids, disease_ids, ortholog_ids)

main_components <- build_radar_components(main_plot_data, "PlotGroup", main_group_order, feature_order, feature_labels)
innocenti_components <- build_radar_components(innocenti_plot_data, "PlotGroup", candidate_group_order, feature_order, feature_labels)
ruzicka_components <- build_radar_components(ruzicka_plot_data, "PlotGroup", candidate_group_order, feature_order, feature_labels)
neal_components <- build_radar_components(neal_plot_data, "PlotGroup", candidate_group_order, feature_order, feature_labels)
overlay_components <- build_overlay_radar_components(
  resolution_membership,
  group_ids = list(
    "SSDOs" = ssdo_ids,
    "Innocenti and Morrow SA candidates" = innocenti_candidate_ids,
    "Ruzicka et al. SA candidates" = ruzicka_candidate_ids,
    "Neal, Truter et al. SA candidates" = neal_candidate_ids
  ),
  group_order = overlay_group_order,
  feature_order = feature_order,
  feature_labels = feature_labels
)

png(
  filename = sup_figure_path("FigS6_radar_context.png"),
  width = 3600,
  height = 3900,
  res = 300,
  bg = "white"
)

layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(15, 15, 2.3))
par(family = "sans", xpd = NA)

par(mar = c(2.8, 2.8, 3.3, 2.8))
draw_radar_panel(main_components, "SSDOs (Figure 4)", "A", main_palette)

par(mar = c(2.8, 2.8, 3.3, 2.8))
draw_radar_panel(
  innocenti_components,
  expression("Innocenti and Morrow (2010)"),
  "B",
  candidate_palette
)

par(mar = c(2.8, 2.8, 3.3, 2.8))
draw_radar_panel(
  ruzicka_components,
  expression("Ruzicka " * italic("et al.") * " (2019)"),
  "C",
  candidate_palette
)

par(mar = c(2.8, 2.8, 3.3, 2.8))
draw_radar_panel(
  neal_components,
  expression("Neal, Truter, " * italic("et al.") * " (unpublished)"),
  "D",
  candidate_palette
)

par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = c("SSDOs / Candidate SA genes", "Disease orthologues", "Orthologues", "Background"),
  horiz = FALSE,
  ncol = 2,
  bty = "n",
  lty = 1,
  lwd = 5,
  col = c("#D55E00", "#0072B2", "#009E73", "#999999"),
  cex = 1.2,
  seg.len = 3,
  x.intersp = 0.9,
  y.intersp = 1.2
)

dev.off()

png(
  filename = sup_figure_path("FigS7_SSDO_candidate_overlay.png"),
  width = 2700,
  height = 2400,
  res = 300,
  bg = "white"
)

par(family = "sans", mar = c(4.6, 2.8, 3.2, 2.8), xpd = NA)

draw_radar_panel(
  overlay_components,
  "",
  "",
  overlay_palette
)

legend(
  "bottom",
  legend = overlay_group_order,
  horiz = FALSE,
  bty = "n",
  ncol = 2,
  lty = 1,
  lwd = 5,
  col = overlay_palette[overlay_group_order],
  cex = 1.15,
  seg.len = 3,
  x.intersp = 0.8,
  y.intersp = 1.1,
  inset = c(0, -0.13)
)

dev.off()
