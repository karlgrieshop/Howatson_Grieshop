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
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
figure_path <- function(...) file.path(project_root, "Figures", ...)


suppressPackageStartupMessages({
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
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = base_size + 2, margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}


# Figure 1 is exported at 12 x 10 inches as a two-column figure.
# These single-column stacked exports are sized so that, after manuscript
# scaling, their text remains close to the apparent size used in Figure 1.
single_column_width <- 6.8
two_panel_height <- 10
three_panel_height <- 14.5


read_summary_data <- function(filename) {
  read.csv(summary_path(filename), stringsAsFactors = FALSE) |>
    transform(GeneType = factor(GeneType, levels = ordered_gene_types))
}


plot_metric_panel <- function(summary_data, y_label, tag, show_x_labels = FALSE) {
  p <- ggplot(summary_data, aes(x = GeneType, y = q50, colour = GeneType)) +
    geom_pointrange(aes(ymin = q05, ymax = q95), linewidth = 1.2) +
    geom_text(aes(y = q95, label = SampleSize), size = 3.9, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = gene_type_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = NULL, y = y_label, tag = tag) +
    publication_theme()

  if (show_x_labels) {
    p + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  } else {
    p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
}


rmf_data <- read_summary_data("rMF Plot Data.csv")
tissue_data <- read_summary_data("Tissue Specificity Plot Data.csv")
rlnd_data <- read_summary_data("rlnd Plot Data.csv")
tajd_data <- read_summary_data("tajD Plot Data.csv")
dos_data <- read_summary_data("DoS Plot Data.csv")


figure2 <- plot_metric_panel(
  summary_data = rmf_data,
  y_label = expression(italic(r)[italic(MF)]),
  tag = "A",
  show_x_labels = FALSE
) /
  plot_metric_panel(
    summary_data = tissue_data,
    y_label = expression("Tissue specificity (" * tau * ")"),
    tag = "B",
    show_x_labels = TRUE
  )


figure3 <- plot_metric_panel(
  summary_data = rlnd_data,
  y_label = expression(tilde(pi)[italic(r) * italic(N)]),
  tag = "A",
  show_x_labels = FALSE
) /
  plot_metric_panel(
    summary_data = dos_data,
    y_label = expression(italic(DoS)),
    tag = "B",
    show_x_labels = FALSE
  ) /
  plot_metric_panel(
    summary_data = tajd_data,
    y_label = expression("Tajima's " * italic(D)),
    tag = "C",
    show_x_labels = TRUE
  )


ggsave(
  filename = figure_path("Fig2_rmf_tau.png"),
  plot = figure2,
  width = single_column_width,
  height = two_panel_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = figure_path("Fig3_rlns_tajD_DoS.png"),
  plot = figure3,
  width = single_column_width,
  height = three_panel_height,
  dpi = 300,
  bg = "white"
)
