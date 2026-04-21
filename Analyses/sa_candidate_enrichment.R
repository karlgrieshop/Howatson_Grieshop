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
        "Could not locate the project root. Run this script from inside the Final_Year_Project_github repository.",
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

  data.frame(Symbol = symbols, stringsAsFactors = FALSE) |>
    distinct() |>
    left_join(symbol_map, by = c("Symbol" = "Fly.Symbol"))
}


load_innocenti_candidates <- function(gene_reference) {
  innocenti_raw <- read.csv(
    data_path("unused", "legacy_data", "Innocenti_and_Morrow_ANTAGONISTIC.csv"),
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


load_ruzicka_candidates <- function() {
  read.csv(
    data_path("unused", "legacy_data", "Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct()
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


run_enrichment_test <- function(universe_ids, ssdo_ids, candidate_ids, study_name, source_note) {
  analysis_df <- data.frame(FlyBaseID = universe_ids, stringsAsFactors = FALSE) |>
    mutate(
      IsSSDO = FlyBaseID %in% ssdo_ids,
      IsCandidate = FlyBaseID %in% candidate_ids
    )

  contingency <- table(analysis_df$IsSSDO, analysis_df$IsCandidate)
  fisher_result <- fisher.test(contingency)

  counts <- as.data.frame.matrix(contingency)
  counts$IsSSDO <- rownames(counts)
  rownames(counts) <- NULL

  summary_row <- data.frame(
    Study = study_name,
    Source = source_note,
    UniverseSize = length(universe_ids),
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
      Label = paste0(Candidate, "/", Candidate + NonCandidate)
    )

  list(summary = summary_row, plot_data = plot_df, overlap = analysis_df |> filter(IsSSDO, IsCandidate))
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

innocenti_result <- run_enrichment_test(
  universe_ids = protein_coding_universe,
  ssdo_ids = ssdo_ids,
  candidate_ids = innocenti$ids$FlyBaseID,
  study_name = "Innocenti and Morrow",
  source_note = "Mapped from legacy antagonistic-gene symbols"
)

ruzicka_result <- run_enrichment_test(
  universe_ids = protein_coding_universe,
  ssdo_ids = ssdo_ids,
  candidate_ids = ruzicka_ids$FlyBaseID,
  study_name = "Ruzicka et al.",
  source_note = "Legacy Ruzicka candidate gene table"
)

neal_result <- run_enrichment_test(
  universe_ids = protein_coding_universe,
  ssdo_ids = ssdo_ids,
  candidate_ids = neal_ids$FlyBaseID,
  study_name = "Neal_Truter_Grieshop",
  source_note = "gene_hits_snps_90.csv"
)

summary_table <- bind_rows(
  innocenti_result$summary,
  ruzicka_result$summary,
  neal_result$summary
) |>
  mutate(PAdjustBH = p.adjust(PValue, method = "BH"))

plot_data <- bind_rows(
  innocenti_result$plot_data,
  ruzicka_result$plot_data,
  neal_result$plot_data
)

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
  bind_rows(
    mutate(innocenti_result$overlap, Study = "Innocenti and Morrow"),
    mutate(ruzicka_result$overlap, Study = "Ruzicka et al."),
    mutate(neal_result$overlap, Study = "Neal_Truter_Grieshop")
  ),
  data_path("Summary_tables", "SA Candidate Enrichment Overlap IDs.csv"),
  row.names = FALSE
)

sa_enrichment_plot <- ggplot(plot_data, aes(x = Group, y = ProportionCandidate, fill = Group)) +
  geom_col(width = 0.6, colour = "black") +
  geom_text(aes(label = Label), vjust = -0.5, size = 3.8) +
  facet_wrap(~ Study, nrow = 1) +
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
  width = 10,
  height = 4.5,
  dpi = 300
)
