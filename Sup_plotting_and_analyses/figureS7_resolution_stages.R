#!/usr/bin/env Rscript
# FigS7: resolution-stage rMF and SA effects, with within-category overlap-aware tests.
find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "Data", "outputs", "sa_resolution_feature_membership.csv"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) stop("Run from inside the SSDO repository.", call. = FALSE)
    current <- parent
  }
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(patchwork)
})
project_root <- find_project_root()
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)
input_path <- function(...) file.path(project_root, "Data", "inputs", ...)
output_path <- function(...) file.path(project_root, "Data", "outputs", ...)
gene_order <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
palette <- c("SSDOs" = "#D55E00", "Disease orthologues" = "#0072B2",
             "Orthologues" = "#009E73", "Background" = "#999999")
stage_order <- c("Penultimate", "Ultimate")
membership_flags <- c("MB", "TS", "HMB", "HTS")
gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_ids <- unique(gene_reference$FlyBaseID[
  !is.na(gene_reference$FlyBaseID) & grepl("^FBgn", gene_reference$FlyBaseID) &
    gene_reference$gene_biotype %in% "protein_coding"
])
normalise_chromosome <- function(x) trimws(gsub("^chr", "", as.character(x), ignore.case = TRUE))
genes <- gene_reference |>
  filter(gene_biotype == "protein_coding", !is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
  transmute(FlyBaseID, chromosome = normalise_chromosome(chromosome),
            start = as.integer(start), end = as.integer(end)) |>
  filter(!is.na(chromosome), nzchar(chromosome), !is.na(start), !is.na(end)) |>
  distinct()
snps <- fread(input_path("Assoc_Data.txt"), select = c("Chromosome", "Basepair", "Wald_Stat")) |>
  as_tibble() |>
  transmute(chromosome = normalise_chromosome(Chromosome), start = as.integer(Basepair),
            end = as.integer(Basepair), Value = as.numeric(Wald_Stat)) |>
  filter(!is.na(chromosome), nzchar(chromosome), !is.na(start), is.finite(Value))
genes_dt <- as.data.table(genes)
snps_dt <- as.data.table(snps)
setkey(genes_dt, chromosome, start, end)
setkey(snps_dt, chromosome, start, end)
ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"))$FlyBaseID
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"))$FlyBaseID
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"))$FlyBaseID
assign_gene_types <- function(values) {
  values |>
    mutate(GeneType = case_when(
      FlyBaseID %in% ssdo_ids ~ "SSDOs",
      FlyBaseID %in% disease_ids ~ "Disease orthologues",
      FlyBaseID %in% ortholog_ids ~ "Orthologues",
      TRUE ~ "Background"
    ))
}
gene_values <- as_tibble(foverlaps(snps_dt, genes_dt, nomatch = 0L)) |>
  group_by(FlyBaseID) |>
  summarise(Value = max(Value), .groups = "drop") |>
  mutate(Source = "Ruzicka et al. from Assoc_Data.txt (max Wald_Stat)") |>
  assign_gene_types()
rmf_source <- read.csv(input_path("SBGEandSSSdataForMBE.csv"), stringsAsFactors = FALSE)
names(rmf_source)[2] <- "FlyBaseID"
rmf_values <- rmf_source |>
  distinct(FlyBaseID, .keep_all = TRUE) |>
  transmute(FlyBaseID, Value = as.numeric(rmf),
            Source = "SBGEandSSSdataForMBE.csv (rmf)") |>
  filter(FlyBaseID %in% protein_coding_ids, is.finite(Value)) |>
  assign_gene_types()
membership <- read.csv(
  file.path(project_root, "Data", "outputs", "sa_resolution_feature_membership.csv"),
  stringsAsFactors = FALSE
)
if (!all(c("FlyBaseID", membership_flags) %in% names(membership)) ||
    anyDuplicated(membership$FlyBaseID) ||
    !all(vapply(membership[membership_flags], is.logical, logical(1)))) {
  stop("SA-resolution table must have unique gene IDs and logical MB/TS/HMB/HTS flags.", call. = FALSE)
}
stage_permutation_tests <- function(values, metric, input, n_perm = 10000L, seed = 1L) {
  if (!setequal(unique(values$GeneType), gene_order) ||
      !setequal(unique(values$ResolutionStage), stage_order) ||
      anyDuplicated(values[c("FlyBaseID", "ResolutionStage")]) ||
      anyDuplicated(unique(values[c("FlyBaseID", "GeneType")])$FlyBaseID) ||
      any(!is.finite(values$Value))) {
    stop("Expected four mutually exclusive gene categories and unique finite estimates per stage.", call. = FALSE)
  }
  tests <- bind_rows(lapply(seq_along(gene_order), function(i) {
    category_values <- filter(values, GeneType == gene_order[i])
    gene_values <- distinct(category_values, FlyBaseID, Value)
    if (anyDuplicated(gene_values$FlyBaseID)) {
      stop("Shared genes must have the same estimate in both stages.", call. = FALSE)
    }
    value_map <- setNames(gene_values$Value, gene_values$FlyBaseID)
    group1_ids <- category_values$FlyBaseID[category_values$ResolutionStage == "Penultimate"]
    group2_ids <- category_values$FlyBaseID[category_values$ResolutionStage == "Ultimate"]
    if (length(group1_ids) < 2L || length(group2_ids) < 2L) {
      stop("Insufficient genes for ", metric, " / ", gene_order[i], call. = FALSE)
    }
    shared_ids <- intersect(group1_ids, group2_ids)
    group1_only <- setdiff(group1_ids, shared_ids)
    group2_only <- setdiff(group2_ids, shared_ids)
    shared_values <- unname(value_map[shared_ids])
    pooled_exclusive <- unname(value_map[c(group1_only, group2_only)])
    n1 <- length(group1_only)
    n2 <- length(group2_only)
    mean1 <- mean(value_map[group1_ids])
    mean2 <- mean(value_map[group2_ids])
    observed <- mean1 - mean2
    category_seed <- seed + i - 1L
    set.seed(category_seed)
    # Table S2's null: shared genes stay fixed in both groups; shuffle exclusive members.
    permuted <- replicate(n_perm, {
      shuffled <- if (length(pooled_exclusive)) {
        pooled_exclusive[sample.int(length(pooled_exclusive))]
      } else numeric(0)
      mean(c(shared_values, shuffled[seq_len(n1)])) -
        mean(c(shared_values, shuffled[n1 + seq_len(n2)]))
    })
    data.frame(
      Metric = metric, GeneGroup = gene_order[i], Group1 = "Penultimate", Group2 = "Ultimate",
      Group1Definition = "MB or TS", Group2Definition = "HMB or HTS",
      NGroup1 = length(group1_ids), NGroup2 = length(group2_ids),
      NOverlap = length(shared_ids), NGroup1Only = n1, NGroup2Only = n2,
      NUniqueGenes = length(value_map),
      MeanGroup1 = mean1, MeanGroup2 = mean2, MeanDifference = observed,
      PValue = (sum(abs(permuted) >= abs(observed)) + 1) / (n_perm + 1),
      Permutations = n_perm, Seed = category_seed, Alternative = "two-sided",
      Method = "Shared genes fixed in both groups; exclusive gene labels permuted at fixed group sizes",
      Universe = "Within each mutually exclusive gene category: protein-coding genes with finite metric values and positive MB/TS or HMB/HTS membership; other genes excluded",
      Input = input, stringsAsFactors = FALSE
    )
  }))
  tests$PAdjustedBH <- p.adjust(tests$PValue, method = "BH")
  tests$AdjustmentFamily <- paste(metric, ": four within-category two-sided Penultimate-versus-Ultimate tests")
  tests
}
build_panel <- function(metric_values, y_label, prefix, values_source) {
  if (nrow(metric_values) == 0L || anyDuplicated(metric_values$FlyBaseID) ||
      any(!is.finite(metric_values$Value))) {
    stop("Expected unique finite gene estimates for ", prefix, call. = FALSE)
  }
  audit <- metric_values |>
    left_join(membership |> select(FlyBaseID, all_of(membership_flags)), by = "FlyBaseID") |>
    mutate(
      Penultimate = MB %in% TRUE | TS %in% TRUE,
      Ultimate = HMB %in% TRUE | HTS %in% TRUE,
      IncludedInPlot = Penultimate | Ultimate,
      ValuesSource = values_source,
      CategorySource = "Data/outputs/sexspecific_disease_ids.csv; all_disease_ortholog_ids.csv; all_ortholog_ids.csv (priority order)",
      MembershipSource = "Data/outputs/sa_resolution_feature_membership.csv"
    )
  plot_values <- bind_rows(lapply(stage_order, function(stage) {
    audit |>
      filter(.data[[stage]]) |>
      transmute(FlyBaseID, ResolutionStage = stage, GeneType, Value, Source)
  })) |>
    mutate(ResolutionStage = factor(ResolutionStage, levels = stage_order),
           GeneType = factor(GeneType, levels = gene_order)) |>
    arrange(ResolutionStage, GeneType, FlyBaseID)
  set.seed(1L)
  summary_data <- bind_rows(lapply(stage_order, function(stage) {
    bind_rows(lapply(gene_order, function(gene_type) {
      x <- plot_values$Value[plot_values$ResolutionStage == stage & plot_values$GeneType == gene_type]
      if (length(x) < 2L) stop("Insufficient genes in ", stage, " / ", gene_type, call. = FALSE)
      boot <- replicate(5000, mean(sample(x, size = length(x), replace = TRUE)))
      data.frame(
        ResolutionStage = stage, GeneType = gene_type, Mean = mean(x),
        q05 = unname(quantile(boot, 0.025)), q50 = unname(quantile(boot, 0.5)),
        q95 = unname(quantile(boot, 0.975)), SampleSize = length(x),
        BootstrapReplicates = 5000L, Seed = 1L
      )
    }))
  })) |>
    mutate(ResolutionStage = factor(ResolutionStage, levels = stage_order),
           GeneType = factor(GeneType, levels = gene_order))
  write.csv(audit, summary_path(paste0(prefix, "gene_inclusion_audit.csv")), row.names = FALSE)
  write.csv(plot_values, summary_path(paste0(prefix, "gene_level_values.csv")), row.names = FALSE)
  write.csv(summary_data, summary_path(paste0(prefix, "bootstrap_summaries.csv")), row.names = FALSE)

  dodge <- position_dodge(width = 0.65)
  figure <- ggplot(summary_data, aes(ResolutionStage, q50, colour = GeneType, group = GeneType)) +
    geom_hline(yintercept = mean(metric_values$Value), colour = "#999999", linewidth = 0.6) +
    geom_pointrange(aes(ymin = q05, ymax = q95), position = dodge, linewidth = 1) +
    geom_text(aes(y = q95, label = SampleSize),
              position = dodge, size = 3.5, colour = "black", vjust = -0.75) +
    scale_colour_manual(values = palette, breaks = gene_order, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    scale_x_discrete(labels = c(Penultimate = "Penultimate\n(MB or TS)", Ultimate = "Ultimate\n(HMB or HTS)")) +
    labs(x = NULL, y = y_label, colour = NULL) +
    guides(colour = guide_legend(ncol = 2)) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 16),
      plot.tag = element_text(size = 18, face = "bold"),
      axis.line = element_blank(), legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
  list(figure = figure, audit = audit, plot_values = plot_values, summary_data = summary_data)
}
sa_panel <- build_panel(
  gene_values, "SA effect", "FigS7_sa_effect_",
  "Data/inputs/Assoc_Data.txt; Data/outputs/dmel_gene_reference.csv (protein-coding)"
)
rmf_panel <- build_panel(
  rmf_values, expression(italic(r)[MF]), "FigS7_rmf_",
  "Data/inputs/SBGEandSSSdataForMBE.csv (rmf); Data/outputs/dmel_gene_reference.csv (protein-coding)"
)
reference_means <- bind_rows(lapply(c("rMF", "SA effect"), function(metric) {
  values <- if (metric == "rMF") rmf_values else gene_values
  data.frame(
    Metric = metric, Mean = mean(values$Value), NGenes = nrow(values),
    Universe = "All protein-coding genes with finite metric values; all gene categories; no resolution-stage restriction",
    Source = unique(values$Source), stringsAsFactors = FALSE
  )
}))
write.csv(reference_means, summary_path("FigS7_genome_wide_reference_means.csv"), row.names = FALSE)
sa_tests <- stage_permutation_tests(
  sa_panel$plot_values, "SA effect", "Data/Summary_tables/FigS7_sa_effect_gene_level_values.csv"
)
rmf_tests <- stage_permutation_tests(
  rmf_panel$plot_values, "rMF", "Data/Summary_tables/FigS7_rmf_gene_level_values.csv"
)
write.csv(sa_tests, summary_path("FigS7_sa_effect_permutation_tests.csv"), row.names = FALSE)
write.csv(rmf_tests, summary_path("FigS7_rmf_permutation_tests.csv"), row.names = FALSE)
figure <- (rmf_panel$figure + sa_panel$figure) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
ggsave(sup_figure_path("FigS7_resolution_stages.png"), figure,
       width = 12, height = 5.5, dpi = 300, bg = "white")
message("Plotted rMF (A) and SA effects (B), with metric-specific genome-wide means.")
