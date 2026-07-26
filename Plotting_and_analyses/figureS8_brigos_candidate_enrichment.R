# Brigos candidate-gene enrichment analyses (Sup Fig S1 style).
# Uses Brigos-Baril strict FlyBase IDs as the candidate SA set.
# Keeps the same study-matched universes logic as sa_candidate_enrichment.R:
#   - Innocenti & Morrow (Affx assay universe)
#   - Ruzicka et al. (Assoc_Data SNP-mapped universe)
#   - Neal (2024) (true_cmh_results SNP-mapped universe)
#
# Runs focal tests for:
#   1) SSDOs
#   2) Disease orthologues
#   3) Orthologues
#
# Outputs:
#   - Summary table with OR, P, CI, BH-adjusted P
#   - Overlap table (FBgns)
#   - S1-style barplot

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
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
})

dir.create(summary_path(), recursive = TRUE, showWarnings = FALSE)
dir.create(sup_figure_path(), recursive = TRUE, showWarnings = FALSE)

simple_theme <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1, face = "plain"),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "none",
      axis.text = element_text(size = base_size, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      strip.text.x = ggtext::element_markdown(size = base_size + 1),
      plot.margin = margin(10, 10, 10, 10)
    )
}

map_symbols_to_flybase <- function(symbols, gene_reference) {
  symbol_map <- gene_reference |>
    dplyr::filter(!is.na(Fly.Symbol), nzchar(Fly.Symbol)) |>
    dplyr::select(FlyBaseID, Fly.Symbol) |>
    dplyr::distinct() |>
    dplyr::add_count(Fly.Symbol, name = "n_ids") |>
    dplyr::filter(n_ids == 1) |>
    dplyr::select(-n_ids)
  
  data.frame(Symbol = as.character(symbols), stringsAsFactors = FALSE) |>
    dplyr::distinct() |>
    dplyr::left_join(symbol_map, by = c("Symbol" = "Fly.Symbol"))
}

normalise_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  trimws(x)
}

map_snp_positions_to_genes <- function(snp_positions, gene_reference) {
  genes <- gene_reference |>
    dplyr::transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(start),
      end = as.integer(end),
      FlyBaseID = as.character(FlyBaseID),
      gene_biotype = as.character(gene_biotype)
    ) |>
    dplyr::filter(
      !is.na(FlyBaseID),
      grepl("^FBgn", FlyBaseID),
      !is.na(chromosome),
      nzchar(chromosome),
      !is.na(start),
      !is.na(end)
    ) |>
    dplyr::distinct(FlyBaseID, chromosome, start, end, gene_biotype)
  
  snps <- snp_positions |>
    dplyr::transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(position),
      end = as.integer(position)
    ) |>
    dplyr::filter(!is.na(chromosome), nzchar(chromosome), !is.na(start)) |>
    dplyr::distinct()
  
  genes_dt <- data.table::as.data.table(genes)
  snps_dt <- data.table::as.data.table(snps)
  
  data.table::setkey(genes_dt, chromosome, start, end)
  data.table::setkey(snps_dt, chromosome, start, end)
  
  overlaps <- data.table::foverlaps(snps_dt, genes_dt, nomatch = 0L)
  
  tibble::as_tibble(overlaps) |>
    dplyr::transmute(FlyBaseID = as.character(FlyBaseID)) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::distinct()
}

load_innocenti_candidates <- function(gene_reference) {
  innocenti_raw <- read.csv(
    input_path("Innocenti_and_Morrow_ANTAGONISTIC.csv"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::transmute(Symbol = as.character(Symbol)) |>
    dplyr::filter(!is.na(Symbol), nzchar(Symbol))
  
  mapped <- map_symbols_to_flybase(unique(innocenti_raw$Symbol), gene_reference)
  
  list(
    ids = mapped |>
      dplyr::filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
      dplyr::transmute(FlyBaseID) |>
      dplyr::distinct(),
    mapping = mapped
  )
}

load_innocenti_matched_universe <- function(gene_reference) {
  affx <- read.csv(
    output_path("Affx_Dros_2.0.csv"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::transmute(Symbol = as.character(Symbol)) |>
    dplyr::filter(!is.na(Symbol), nzchar(Symbol)) |>
    dplyr::distinct()
  
  map_symbols_to_flybase(affx$Symbol, gene_reference) |>
    dplyr::filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
    dplyr::transmute(FlyBaseID) |>
    dplyr::distinct()
}

load_ruzicka_candidates <- function() {
  read.csv(
    input_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    dplyr::transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::distinct()
}

load_ruzicka_matched_universe <- function(gene_reference) {
  assoc <- read.delim(
    input_path("Assoc_Data.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    dplyr::transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair)
    )
  
  map_snp_positions_to_genes(assoc, gene_reference)
}

load_neal_candidates <- function() {
  read.csv(
    input_path("gene_hits_snps_90.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    dplyr::transmute(FlyBaseID = as.character(ensembl_gene_id)) |>
    dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    dplyr::distinct()
}

load_neal_matched_universe <- function(gene_reference) {
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
    dplyr::transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair)
    )
  
  map_snp_positions_to_genes(snps, gene_reference)
}

run_enrichment_test <- function(
    universe_ids,
    focal_ids,
    candidate_ids,
    study_name,
    universe_type,
    source_note,
    focal_group,
    panel_label
) {
  analysis_df <- data.frame(FlyBaseID = unique(universe_ids), stringsAsFactors = FALSE) |>
    dplyr::mutate(
      IsFocal = FlyBaseID %in% focal_ids,
      IsCandidate = FlyBaseID %in% candidate_ids
    )
  
  contingency <- table(analysis_df$IsFocal, analysis_df$IsCandidate)
  fisher_result <- fisher.test(contingency)
  
  summary_row <- data.frame(
    FocalGroup = focal_group,
    Study = study_name,
    UniverseType = universe_type,
    Source = source_note,
    UniverseSize = nrow(analysis_df),
    FocalCount = sum(analysis_df$IsFocal),
    CandidateCount = sum(analysis_df$IsCandidate),
    OverlapCount = sum(analysis_df$IsFocal & analysis_df$IsCandidate),
    OddsRatio = unname(fisher_result$estimate),
    PValue = fisher_result$p.value,
    ConfLow = fisher_result$conf.int[1],
    ConfHigh = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )
  
  plot_df <- analysis_df |>
    dplyr::mutate(
      DisplayGroup = ifelse(IsFocal, "Focal genes", "Other genes"),
      FillGroup = ifelse(IsFocal, focal_group, "Other genes")
    ) |>
    dplyr::group_by(DisplayGroup, FillGroup) |>
    dplyr::summarise(
      Candidate = sum(IsCandidate),
      NonCandidate = sum(!IsCandidate),
      ProportionCandidate = Candidate / (Candidate + NonCandidate),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      FocalGroup = focal_group,
      PanelLabel = panel_label,
      Study = study_name,
      UniverseType = universe_type,
      Label = paste0(Candidate, "/", Candidate + NonCandidate)
    )
  
  overlap <- analysis_df |>
    dplyr::filter(IsFocal, IsCandidate) |>
    dplyr::mutate(FocalGroup = focal_group, Study = study_name, UniverseType = universe_type)
  
  list(summary = summary_row, plot_data = plot_df, overlap = overlap)
}

# ---- Main execution block -------------------------------------------------

gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  dplyr::pull(FlyBaseID) |>
  unique()

disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  dplyr::pull(FlyBaseID) |>
  unique()

ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  dplyr::pull(FlyBaseID) |>
  unique()

# Brigos candidate set (strict)
brigos_ids <- read.csv(output_path("brigos_barril_ids_strict.csv"), stringsAsFactors = FALSE) |>
  dplyr::transmute(FlyBaseID = as.character(FlyBaseID)) |>
  dplyr::filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
  dplyr::distinct() |>
  dplyr::pull(FlyBaseID)

disease_only_ids <- setdiff(disease_ids, ssdo_ids)
ortholog_only_ids <- setdiff(ortholog_ids, disease_ids)

study_order <- c("Innocenti and Morrow", "Ruzicka et al.", "Neal_2024")
study_display_order <- c(
  "Innocenti and Morrow (2010)",
  "Ruzicka <i>et al.</i> (2019)",
  "Neal (2024)"
)
study_display_map <- c(
  "Innocenti and Morrow" = "Innocenti and Morrow (2010)",
  "Ruzicka et al." = "Ruzicka <i>et al.</i> (2019)",
  "Neal_2024" = "Neal (2024)"
)

focal_order <- c("SSDOs", "Disease orthologues", "Orthologues")
panel_order <- c(
  "Brigos genes among SSDOs",
  "Brigos genes among Disease Orthologues",
  "Brigos genes among Orthologues"
)

# Study-matched universes are unchanged
matched_universes <- list(
  "Innocenti and Morrow" = load_innocenti_matched_universe(gene_reference)$FlyBaseID,
  "Ruzicka et al." = load_ruzicka_matched_universe(gene_reference)$FlyBaseID,
  "Neal_2024" = load_neal_matched_universe(gene_reference)$FlyBaseID
)

# Candidate sets per study remain the original SA-candidate definitions
innocenti_candidates <- load_innocenti_candidates(gene_reference)$ids$FlyBaseID
ruzicka_candidates <- load_ruzicka_candidates()$FlyBaseID
neal_candidates <- load_neal_candidates()$FlyBaseID

# Intersect Brigos IDs with study universes to keep assay-matched framing
brigos_in_study <- list(
  "Innocenti and Morrow" = intersect(brigos_ids, matched_universes[["Innocenti and Morrow"]]),
  "Ruzicka et al." = intersect(brigos_ids, matched_universes[["Ruzicka et al."]]),
  "Neal_2024" = intersect(brigos_ids, matched_universes[["Neal_2024"]])
)

# ---- Main analysis workflow -----------------------------------------------

ssdo_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = ssdo_ids,
    candidate_ids = brigos_in_study[["Innocenti and Morrow"]],
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "SSDOs",
    panel_label = "Brigos genes among SSDOs"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = ssdo_ids,
    candidate_ids = brigos_in_study[["Ruzicka et al."]],
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "SSDOs",
    panel_label = "Brigos genes among SSDOs"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = ssdo_ids,
    candidate_ids = brigos_in_study[["Neal_2024"]],
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "SSDOs",
    panel_label = "Brigos genes among SSDOs"
  )
)

disease_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = disease_only_ids,
    candidate_ids = brigos_in_study[["Innocenti and Morrow"]],
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "Disease orthologues",
    panel_label = "Brigos genes among Disease Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = disease_only_ids,
    candidate_ids = brigos_in_study[["Ruzicka et al."]],
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "Disease orthologues",
    panel_label = "Brigos genes among Disease Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = disease_only_ids,
    candidate_ids = brigos_in_study[["Neal_2024"]],
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "Disease orthologues",
    panel_label = "Brigos genes among Disease Orthologues"
  )
)

orthologue_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = ortholog_only_ids,
    candidate_ids = brigos_in_study[["Innocenti and Morrow"]],
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "Orthologues",
    panel_label = "Brigos genes among Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = ortholog_only_ids,
    candidate_ids = brigos_in_study[["Ruzicka et al."]],
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "Orthologues",
    panel_label = "Brigos genes among Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = ortholog_only_ids,
    candidate_ids = brigos_in_study[["Neal_2024"]],
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "Orthologues",
    panel_label = "Brigos genes among Orthologues"
  )
)

all_results <- c(ssdo_results, disease_results, orthologue_results)

summary_table <- dplyr::bind_rows(lapply(all_results, `[[`, "summary")) |>
  dplyr::mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    Study = factor(Study, levels = study_order)
  ) |>
  dplyr::arrange(FocalGroup, Study) |>
  dplyr::mutate(PAdjustBH = p.adjust(PValue, method = "BH"))

plot_data <- dplyr::bind_rows(lapply(all_results, `[[`, "plot_data")) |>
  dplyr::mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    PanelLabel = factor(PanelLabel, levels = panel_order),
    Study = factor(Study, levels = study_order),
    StudyDisplay = factor(dplyr::recode(as.character(Study), !!!study_display_map), levels = study_display_order),
    DisplayGroup = factor(DisplayGroup, levels = c("Focal genes", "Other genes"))
  )

overlap_table <- dplyr::bind_rows(lapply(all_results, `[[`, "overlap")) |>
  dplyr::mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    Study = factor(Study, levels = study_order)
  ) |>
  dplyr::arrange(FocalGroup, Study, FlyBaseID)

write.csv(
  summary_table,
  summary_path("TableS_Brigos_Candidate_Enrichment_Tests.csv"),
  row.names = FALSE
)

write.csv(
  overlap_table,
  summary_path("TableS_Brigos_Candidate_Enrichment_Overlaps.csv"),
  row.names = FALSE
)

brigos_enrichment_plot <- ggplot(plot_data, aes(x = DisplayGroup, y = ProportionCandidate, fill = FillGroup)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = Label), vjust = -0.5, size = 4) +
  geom_text(
    data = summary_table |>
      dplyr::mutate(
        PanelLabel = factor(
          dplyr::case_when(
            FocalGroup == "SSDOs" ~ "Brigos genes among SSDOs",
            FocalGroup == "Disease orthologues" ~ "Brigos genes among Disease Orthologues",
            TRUE ~ "Brigos genes among Orthologues"
          ),
          levels = panel_order
        ),
        StudyDisplay = factor(dplyr::recode(as.character(Study), !!!study_display_map), levels = study_display_order)
      ) |>
      dplyr::filter(PValue < 0.05) |>
      dplyr::mutate(x = 1.5, y = max(plot_data$ProportionCandidate) * 1.12, label = "*"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 10,
    fontface = "bold"
  ) +
  facet_grid(
    PanelLabel ~ StudyDisplay,
    switch = "y",
    labeller = labeller(PanelLabel = function(x) rep("% Brigos genes", length(x)))
  ) +
  scale_fill_manual(values = c(
    "SSDOs" = "#D55E00",
    "Disease orthologues" = "#0072B2",
    "Orthologues" = "#009E73",
    "Other genes" = "#999999"
  )) +
  scale_y_continuous(
    limits = c(0, max(plot_data$ProportionCandidate) * 1.18),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  simple_theme() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "plain", size = 15)
  )

ggsave(
  filename = sup_figure_path("FigS_Brigos_Candidate_Enrichment.png"),
  plot = brigos_enrichment_plot,
  width = 11,
  height = 10.1,
  dpi = 300
)

message("Done.")
message("Summary table: ", summary_path("TableS_Brigos_Candidate_Enrichment_Tests.csv"))
message("Overlap table: ", summary_path("TableS_Brigos_Candidate_Enrichment_Overlaps.csv"))
message("Figure: ", sup_figure_path("FigS_Brigos_Candidate_Enrichment.png"))