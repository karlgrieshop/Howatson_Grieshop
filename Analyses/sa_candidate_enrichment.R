find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Analyses")) &&
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
data_path <- function(...) file.path(project_root, "Data", ...)
figure_path <- function(...) file.path(project_root, "Figures", ...)


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})


simple_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      legend.position = "none",
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
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


normalise_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  trimws(x)
}


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


load_innocenti_candidates <- function(gene_reference) {
  innocenti_raw <- read.csv(
    data_path("Innocenti_and_Morrow_ANTAGONISTIC.csv"),
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
    data_path("Affx_Dros_2.0.csv"),
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
    data_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
}


load_ruzicka_matched_universe <- function(gene_reference) {
  assoc <- read.delim(
    data_path("Assoc_Data.txt"),
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
    data_path("gene_hits_snps_90.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(ensembl_gene_id)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
}


load_neal_matched_universe <- function(gene_reference) {
  snps <- read.delim(
    data_path("true_cmh_results.txt"),
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


run_enrichment_test <- function(universe_ids, ssdo_ids, candidate_ids, study_name, universe_type, source_note) {
  analysis_df <- data.frame(FlyBaseID = unique(universe_ids), stringsAsFactors = FALSE) |>
    mutate(
      IsSSDO = FlyBaseID %in% ssdo_ids,
      IsCandidate = FlyBaseID %in% candidate_ids
    )

  contingency <- table(analysis_df$IsSSDO, analysis_df$IsCandidate)
  fisher_result <- fisher.test(contingency)

  summary_row <- data.frame(
    Study = study_name,
    UniverseType = universe_type,
    Source = source_note,
    UniverseSize = nrow(analysis_df),
    SSDOCount = sum(analysis_df$IsSSDO),
    CandidateCount = sum(analysis_df$IsCandidate),
    OverlapCount = sum(analysis_df$IsSSDO & analysis_df$IsCandidate),
    OddsRatio = unname(fisher_result$estimate),
    PValue = fisher_result$p.value,
    ConfLow = fisher_result$conf.int[1],
    ConfHigh = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )

  plot_df <- analysis_df |>
    mutate(Group = ifelse(IsSSDO, "SSDOs", "Non-SSDOs")) |>
    group_by(Group) |>
    summarise(
      Candidate = sum(IsCandidate),
      NonCandidate = sum(!IsCandidate),
      ProportionCandidate = Candidate / (Candidate + NonCandidate),
      .groups = "drop"
    ) |>
    mutate(
      Study = study_name,
      UniverseType = universe_type,
      Label = paste0(Candidate, "/", Candidate + NonCandidate)
    )

  overlap <- analysis_df |>
    filter(IsSSDO, IsCandidate) |>
    mutate(Study = study_name, UniverseType = universe_type)

  list(summary = summary_row, plot_data = plot_df, overlap = overlap)
}


gene_reference <- read.csv(data_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
ssdo_ids <- read.csv(data_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |> pull(FlyBaseID) |> unique()
protein_coding_universe <- gene_reference |>
  filter(gene_biotype == "protein_coding") |>
  pull(FlyBaseID) |>
  unique()

innocenti <- load_innocenti_candidates(gene_reference)
ruzicka_ids <- load_ruzicka_candidates()
neal_ids <- load_neal_candidates()

matched_universes <- list(
  "Innocenti and Morrow" = load_innocenti_matched_universe(gene_reference)$FlyBaseID,
  "Ruzicka et al." = load_ruzicka_matched_universe(gene_reference)$FlyBaseID,
  "Neal_Truter_Grieshop" = load_neal_matched_universe(gene_reference)$FlyBaseID
)

broad_results <- list(
  run_enrichment_test(
    universe_ids = protein_coding_universe,
    ssdo_ids = ssdo_ids,
    candidate_ids = innocenti$ids$FlyBaseID,
    study_name = "Innocenti and Morrow",
    universe_type = "Whole fly genome",
    source_note = "Whole-genome / protein-coding background"
  ),
  run_enrichment_test(
    universe_ids = protein_coding_universe,
    ssdo_ids = ssdo_ids,
    candidate_ids = ruzicka_ids$FlyBaseID,
    study_name = "Ruzicka et al.",
    universe_type = "Whole fly genome",
    source_note = "Whole-genome / protein-coding background"
  ),
  run_enrichment_test(
    universe_ids = protein_coding_universe,
    ssdo_ids = ssdo_ids,
    candidate_ids = neal_ids$FlyBaseID,
    study_name = "Neal_Truter_Grieshop",
    universe_type = "Whole fly genome",
    source_note = "Whole-genome / protein-coding background"
  )
)

matched_results <- list(
  run_enrichment_test(
    universe_ids = matched_universes[["Innocenti and Morrow"]],
    ssdo_ids = ssdo_ids,
    candidate_ids = innocenti$ids$FlyBaseID,
    study_name = "Innocenti and Morrow",
    universe_type = "Study-matched universe",
    source_note = "Affx_Dros_2.0 array-assayed genes"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Ruzicka et al."]],
    ssdo_ids = ssdo_ids,
    candidate_ids = ruzicka_ids$FlyBaseID,
    study_name = "Ruzicka et al.",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in Assoc_Data.txt"
  ),
  run_enrichment_test(
    universe_ids = matched_universes[["Neal_Truter_Grieshop"]],
    ssdo_ids = ssdo_ids,
    candidate_ids = neal_ids$FlyBaseID,
    study_name = "Neal_Truter_Grieshop",
    universe_type = "Study-matched universe",
    source_note = "Genes represented by SNPs in true_cmh_results.txt"
  )
)

all_results <- c(broad_results, matched_results)

summary_table <- bind_rows(lapply(all_results, `[[`, "summary")) |>
  group_by(UniverseType) |>
  mutate(PAdjustBH = p.adjust(PValue, method = "BH")) |>
  ungroup()

plot_data <- bind_rows(lapply(all_results, `[[`, "plot_data"))
overlap_table <- bind_rows(lapply(all_results, `[[`, "overlap"))

write.csv(
  summary_table,
  data_path("Summary_tables", "SA Candidate Enrichment Tests.csv"),
  row.names = FALSE
)

write.csv(
  innocenti$mapping,
  data_path("Summary_tables", "Innocenti Candidate Symbol Mapping.csv"),
  row.names = FALSE
)

write.csv(
  overlap_table,
  data_path("Summary_tables", "SA Candidate Enrichment Overlap IDs.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Study = names(matched_universes),
    UniverseSize = vapply(matched_universes, length, integer(1)),
    stringsAsFactors = FALSE
  ),
  data_path("Summary_tables", "SA Candidate Matched Universes.csv"),
  row.names = FALSE
)

sa_enrichment_plot <- ggplot(plot_data, aes(x = Group, y = ProportionCandidate, fill = Group)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = Label), vjust = -0.5, size = 3.6) +
  facet_grid(UniverseType ~ Study) +
  scale_fill_manual(values = c("SSDOs" = "#D55E00", "Non-SSDOs" = "#999999")) +
  scale_y_continuous(
    limits = c(0, max(plot_data$ProportionCandidate) * 1.18),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL,
    y = "Candidate SA genes",
    title = NULL
  ) +
  simple_theme()

ggsave(
  filename = figure_path("SA_candidate_enrichment.pdf"),
  plot = sa_enrichment_plot,
  width = 12,
  height = 7,
  dpi = 300
)
