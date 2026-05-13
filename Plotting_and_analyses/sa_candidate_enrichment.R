# Candidate SA-gene enrichment analyses for FigS1 / TableS1.
# For each study, we define a study-specific assay universe and test whether
# the focal gene class is enriched for candidate SA genes relative to that
# background, exporting both summary statistics and exact overlapping gene IDs.

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
  library(ggtext)
})

dir.create(file.path(project_root, "Sup_figures"), recursive = TRUE, showWarnings = FALSE)


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


# Map published symbols to modern FlyBase IDs only when the mapping is unique.
# This avoids inflating overlaps with ambiguous symbol-to-gene matches.

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


normalise_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  trimws(x)
}



# Convert SNP-level studies into gene-level universes by interval overlap with
# the current FlyBase/Ensembl reference coordinates.
map_snp_positions_to_genes <- function(snp_positions, gene_reference) {
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

  snps <- snp_positions |>
    transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(position),
      end = as.integer(position)
    ) |>
    filter(!is.na(chromosome), nzchar(chromosome), !is.na(start)) |>
    distinct()

  genes_dt <- as.data.table(genes)
  snps_dt <- as.data.table(snps)

  setkey(genes_dt, chromosome, start, end)
  setkey(snps_dt, chromosome, start, end)

  overlaps <- foverlaps(snps_dt, genes_dt, nomatch = 0L)

  as_tibble(overlaps) |>
    transmute(FlyBaseID = as.character(FlyBaseID)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
}


# Innocenti and Morrow requires an extra symbol-mapping step because the source
# study is array/probe oriented rather than already keyed to FlyBase IDs.
load_innocenti_candidates <- function(gene_reference) {
  innocenti_raw <- read.csv(
    input_path("Innocenti_and_Morrow_ANTAGONISTIC.csv"),
    stringsAsFactors = FALSE
  ) |>
    transmute(Symbol = as.character(Symbol)) |>
    filter(!is.na(Symbol), nzchar(Symbol))

  mapped <- map_symbols_to_flybase(unique(innocenti_raw$Symbol), gene_reference)

  list(
    ids = mapped |>
      filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
      transmute(FlyBaseID) |>
      distinct(),
    mapping = mapped
  )
}


load_innocenti_matched_universe <- function(gene_reference) {
  affx <- read.csv(
    output_path("Affx_Dros_2.0.csv"),
    stringsAsFactors = FALSE
  ) |>
    transmute(Symbol = as.character(Symbol)) |>
    filter(!is.na(Symbol), nzchar(Symbol)) |>
    distinct()

  map_symbols_to_flybase(affx$Symbol, gene_reference) |>
    filter(!is.na(FlyBaseID), nzchar(FlyBaseID)) |>
    transmute(FlyBaseID) |>
    distinct()
}


load_ruzicka_candidates <- function() {
  read.csv(
    input_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
}


load_ruzicka_matched_universe <- function(gene_reference) {
  assoc <- read.delim(
    input_path("Assoc_Data.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(
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
    transmute(FlyBaseID = as.character(ensembl_gene_id)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
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
    transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair)
    )

  map_snp_positions_to_genes(snps, gene_reference)
}


# Core enrichment helper: compare the focal gene set against the non-focal
# genes present in the same study-specific assay universe.

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
    mutate(
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
    mutate(
      DisplayGroup = ifelse(IsFocal, "Focal genes", "Other genes"),
      FillGroup = ifelse(IsFocal, focal_group, "Other genes")
    ) |>
    group_by(DisplayGroup, FillGroup) |>
    summarise(
      Candidate = sum(IsCandidate),
      NonCandidate = sum(!IsCandidate),
      ProportionCandidate = Candidate / (Candidate + NonCandidate),
      .groups = "drop"
    ) |>
    mutate(
      FocalGroup = focal_group,
      PanelLabel = panel_label,
      Study = study_name,
      UniverseType = universe_type,
      Label = paste0(Candidate, "/", Candidate + NonCandidate)
    )

  overlap <- analysis_df |>
    filter(IsFocal, IsCandidate) |>
    mutate(FocalGroup = focal_group, Study = study_name, UniverseType = universe_type)

  list(summary = summary_row, plot_data = plot_df, overlap = overlap)
}


# ---- Main execution block -------------------------------------------------

gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
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
  "SA genes among SSDOs",
  "SA genes among Disease Orthologues",
  "SA genes among Orthologues"
)

innocenti <- load_innocenti_candidates(gene_reference)
ruzicka_ids <- load_ruzicka_candidates()
neal_ids <- load_neal_candidates()

# ---- Main analysis workflow -----------------------------------------------

matched_universes <- list(
  "Innocenti and Morrow" = load_innocenti_matched_universe(gene_reference)$FlyBaseID,
  "Ruzicka et al." = load_ruzicka_matched_universe(gene_reference)$FlyBaseID,
  "Neal_2024" = load_neal_matched_universe(gene_reference)$FlyBaseID
)

ssdo_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = ssdo_ids,
    candidate_ids = innocenti$ids$FlyBaseID,
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "SSDOs",
    panel_label = "SA genes among SSDOs"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = ssdo_ids,
    candidate_ids = ruzicka_ids$FlyBaseID,
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "SSDOs",
    panel_label = "SA genes among SSDOs"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = ssdo_ids,
    candidate_ids = neal_ids$FlyBaseID,
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "SSDOs",
    panel_label = "SA genes among SSDOs"
  )
)

disease_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = disease_only_ids,
    candidate_ids = innocenti$ids$FlyBaseID,
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "Disease orthologues",
    panel_label = "SA genes among Disease Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = disease_only_ids,
    candidate_ids = ruzicka_ids$FlyBaseID,
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "Disease orthologues",
    panel_label = "SA genes among Disease Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = disease_only_ids,
    candidate_ids = neal_ids$FlyBaseID,
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "Disease orthologues",
    panel_label = "SA genes among Disease Orthologues"
  )
)

orthologue_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    focal_ids = ortholog_only_ids,
    candidate_ids = innocenti$ids$FlyBaseID,
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes",
    focal_group = "Orthologues",
    panel_label = "SA genes among Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    focal_ids = ortholog_only_ids,
    candidate_ids = ruzicka_ids$FlyBaseID,
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt",
    focal_group = "Orthologues",
    panel_label = "SA genes among Orthologues"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_2024"]],
    focal_ids = ortholog_only_ids,
    candidate_ids = neal_ids$FlyBaseID,
    study_name = "Neal_2024",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt",
    focal_group = "Orthologues",
    panel_label = "SA genes among Orthologues"
  )
)

all_results <- c(ssdo_results, disease_results, orthologue_results)

summary_table <- bind_rows(lapply(all_results, `[[`, "summary")) |>
  mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    Study = factor(Study, levels = study_order)
  ) |>
  arrange(FocalGroup, Study)

plot_data <- bind_rows(lapply(all_results, `[[`, "plot_data")) |>
  mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    PanelLabel = factor(PanelLabel, levels = panel_order),
    Study = factor(Study, levels = study_order),
    StudyDisplay = factor(recode(as.character(Study), !!!study_display_map), levels = study_display_order),
    DisplayGroup = factor(DisplayGroup, levels = c("Focal genes", "Other genes"))
  )

overlap_table <- bind_rows(lapply(all_results, `[[`, "overlap")) |>
  mutate(
    FocalGroup = factor(FocalGroup, levels = focal_order),
    Study = factor(Study, levels = study_order)
  ) |>
  arrange(FocalGroup, Study, FlyBaseID)

write.csv(
  summary_table,
  summary_path("TableS1_SA_Candidate_Enrichment_Tests.csv"),
  row.names = FALSE
)

sa_enrichment_plot <- ggplot(plot_data, aes(x = DisplayGroup, y = ProportionCandidate, fill = FillGroup)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = Label), vjust = -0.5, size = 4) +
  geom_text(
    data = summary_table |>
      mutate(
        PanelLabel = factor(
          case_when(
            FocalGroup == "SSDOs" ~ "SA genes among SSDOs",
            FocalGroup == "Disease orthologues" ~ "SA genes among Disease Orthologues",
            TRUE ~ "SA genes among Orthologues"
          ),
          levels = panel_order
        ),
        StudyDisplay = factor(recode(as.character(Study), !!!study_display_map), levels = study_display_order)
      ) |>
      filter(PValue < 0.05) |>
      mutate(x = 1.5, y = max(plot_data$ProportionCandidate) * 1.12, label = "*"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 10,
    fontface = "bold"
  ) +
  facet_grid(
    PanelLabel ~ StudyDisplay,
    switch = "y",
    labeller = labeller(PanelLabel = function(x) rep("% SA genes", length(x)))
  ) +
  scale_fill_manual(values = c("SSDOs" = "#D55E00", "Disease orthologues" = "#0072B2", "Orthologues" = "#009E73", "Other genes" = "#999999")) +
  scale_y_continuous(
    limits = c(0, max(plot_data$ProportionCandidate) * 1.18),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  simple_theme() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "plain", size = 15)
  )

ggsave(
  filename = sup_figure_path("FigS1_SA_Candidate_Enrichment.png"),
  plot = sa_enrichment_plot,
  width = 11,
  height = 10.1,
  dpi = 300
)
