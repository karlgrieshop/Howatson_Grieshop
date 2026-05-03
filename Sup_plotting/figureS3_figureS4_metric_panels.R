# Assemble the SBGE-stratified supplementary metric figures from the exported
# summary tables produced by `Sup_plotting/sbge_metric_plots.R`.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Sup_plotting")) &&
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
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", "SBGE_supplement", ...)
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)


suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})


ordered_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
ordered_sbge_panels <- c(
  "All",
  "Highly female-biased",
  "Female-biased",
  "Unbiased",
  "Male-biased",
  "Highly male-biased"
)
sbge_axis_labels <- c("All", "Highly FB", "Female", "Unbiased", "Male", "Highly MB")

gene_type_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)


publication_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1),
      axis.line = element_blank(),
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      legend.title = element_text(size = base_size + 2, face = "bold"),
      legend.text = element_text(size = base_size + 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = base_size + 1, margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}


# Read one summary table and immediately impose the plotting order used across
# all SBGE-stratified supplementary panels.
read_summary_data <- function(filename) {
  read.csv(summary_path(filename), stringsAsFactors = FALSE) |>
    transform(
      GeneType = factor(GeneType, levels = ordered_gene_types),
      SexBiasPanel = factor(SexBiasPanel, levels = ordered_sbge_panels)
    )
}


plot_metric_panel_by_sbge <- function(summary_data, y_label, tag, show_x_labels = FALSE) {
  dodge_width <- 0.55

  p <- ggplot(summary_data, aes(x = SexBiasPanel, y = q50, colour = GeneType, group = GeneType)) +
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
    scale_x_discrete(labels = sbge_axis_labels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = NULL, y = y_label, tag = tag, colour = "Gene Type") +
    geom_vline(xintercept = 1.5, colour = "black", linewidth = 1) +
    geom_vline(xintercept = 2.5:5.5, colour = "grey60", linewidth = 0.45) +
    publication_theme()

  if (show_x_labels) {
    p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))
  } else {
    p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
}


# Figure sizes are tuned to keep panel text legible after the multi-panel
# layouts are inserted into the supplement.
two_panel_width <- 12
two_panel_height <- 10
three_panel_width <- 12
three_panel_height <- 14.5


rmf_data <- read_summary_data("rMF Plot Data.csv")
tissue_data <- read_summary_data("Tissue Specificity Plot Data.csv")
rlnd_data <- read_summary_data("rlnd Plot Data.csv")
tajd_data <- read_summary_data("tajD Plot Data.csv")
dos_data <- read_summary_data("DoS Plot Data.csv")


figure_s3 <- (
  plot_metric_panel_by_sbge(
    summary_data = rmf_data,
    y_label = expression(italic(r)[italic(MF)]),
    tag = "A",
    show_x_labels = FALSE
  ) /
    plot_metric_panel_by_sbge(
      summary_data = tissue_data,
      y_label = expression("Tissue specificity (" * tau * ")"),
      tag = "B",
      show_x_labels = TRUE
    )
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")


figure_s4 <- (
  plot_metric_panel_by_sbge(
    summary_data = rlnd_data,
    y_label = expression(tilde(pi)[italic(r) * italic(N)]),
    tag = "A",
    show_x_labels = FALSE
  ) /
    plot_metric_panel_by_sbge(
      summary_data = dos_data,
      y_label = expression(italic(DoS)),
      tag = "B",
      show_x_labels = FALSE
    ) /
    plot_metric_panel_by_sbge(
      summary_data = tajd_data,
      y_label = expression("Tajima's " * italic(D)),
      tag = "C",
      show_x_labels = TRUE
    )
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")


ggsave(
  filename = sup_figure_path("FigS3_rmf_tau.png"),
  plot = figure_s3,
  width = two_panel_width,
  height = two_panel_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = sup_figure_path("FigS4_rlns_tajD_DoS.png"),
  plot = figure_s4,
  width = three_panel_width,
  height = three_panel_height,
  dpi = 300,
  bg = "white"
)
