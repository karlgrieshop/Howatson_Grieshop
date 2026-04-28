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
figure_path <- function(...) file.path(project_root, "Figures", ...)
sup_figure_path <- function(...) file.path(project_root, "Sup_figures", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)


suppressPackageStartupMessages({
  library(dplyr)
})


ordered_gene_types <- c("SSDOs", "Disease orthologues", "Orthologues", "Background")
gene_type_palette <- c(
  "SSDOs" = "#D55E00",
  "Disease orthologues" = "#0072B2",
  "Orthologues" = "#009E73",
  "Background" = "#999999"
)


# Run `Data_wrangling/build_sa_resolution_inputs.R` first. That script assembles
# the genome-wide feature-membership table and records which requested
# SA-resolution features are ready versus still blocked by missing external
# sources.

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
disease_ids <- read.csv(output_path("all_disease_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
ortholog_ids <- read.csv(output_path("all_ortholog_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
resolution_membership_path <- output_path("sa_resolution_feature_membership.csv")
resolution_notes_path <- summary_path("sa_resolution_feature_notes.csv")

if (!file.exists(resolution_membership_path) || !file.exists(resolution_notes_path)) {
  stop(
    paste(
      "Missing SA-resolution input tables.",
      "Run Data_wrangling/build_sa_resolution_inputs.R before running this analysis.",
      sep = " "
    ),
    call. = FALSE
  )
}

resolution_membership <- read.csv(resolution_membership_path, stringsAsFactors = FALSE)
resolution_notes <- read.csv(resolution_notes_path, stringsAsFactors = FALSE)

available_resolution_features <- resolution_notes |>
  filter(Status == "ready") |>
  pull(Feature)
reported_resolution_features <- available_resolution_features

pending_resolution_features <- resolution_notes |>
  filter(Status != "ready") |>
  pull(Feature)


get_resolution_candidate_ids <- function(feature_name) {
  if (!(feature_name %in% available_resolution_features)) {
    stop(
      sprintf(
        "Resolution feature '%s' is not ready. Pending features: %s",
        feature_name,
        paste(pending_resolution_features, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  resolution_membership |>
    filter(.data[[feature_name]]) |>
    pull(FlyBaseID) |>
    unique()
}


run_resolution_enrichment <- function(universe_ids, candidate_ids, study_name, universe_type, source_note) {
  analysis_df <- data.frame(FlyBaseID = unique(universe_ids), stringsAsFactors = FALSE) |>
    mutate(
      IsSSDO = FlyBaseID %in% ssdo_ids,
      IsCandidate = FlyBaseID %in% candidate_ids
    )

  contingency <- table(analysis_df$IsSSDO, analysis_df$IsCandidate)
  fisher_result <- fisher.test(contingency)

  data.frame(
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
}


run_enrichment_family <- function(universe_ids, universe_type_label, source_note_suffix = NULL) {
  results <- bind_rows(lapply(reported_resolution_features, function(feature_name) {
    source_note <- resolution_notes$Source[match(feature_name, resolution_notes$Feature)]
    if (!is.null(source_note_suffix)) {
      source_note <- paste(source_note, source_note_suffix)
    }

    run_resolution_enrichment(
      universe_ids = universe_ids,
      candidate_ids = get_resolution_candidate_ids(feature_name),
      study_name = feature_name,
      universe_type = universe_type_label,
      source_note = source_note
    )
  })) |>
    mutate(PAdjustBH = p.adjust(PValue, method = "BH")) |>
    arrange(PValue, Study)

  overlaps <- bind_rows(lapply(reported_resolution_features, function(feature_name) {
    candidate_ids <- get_resolution_candidate_ids(feature_name)
    overlap <- sort(intersect(intersect(ssdo_ids, universe_ids), candidate_ids))

    data.frame(
      Study = feature_name,
      UniverseType = universe_type_label,
      FlyBaseID = overlap,
      stringsAsFactors = FALSE
    )
  }))

  list(results = results, overlaps = overlaps)
}


two_sample_permutation_test <- function(x, y, n_perm = 10000, seed = 1L) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(
      data.frame(
        obs_diff = NA_real_,
        n_group1 = length(x),
        n_group2 = length(y),
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  set.seed(seed)

  n_x <- length(x)
  pooled <- c(x, y)
  obs_diff <- mean(x) - mean(y)

  perm_diffs <- replicate(n_perm, {
    shuffled <- sample(pooled, size = length(pooled), replace = FALSE)
    mean(shuffled[seq_len(n_x)]) - mean(shuffled[-seq_len(n_x)])
  })

  data.frame(
    obs_diff = obs_diff,
    n_group1 = length(x),
    n_group2 = length(y),
    p_value = (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1),
    stringsAsFactors = FALSE
  )
}


assign_resolution_groups <- function(data) {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))

  data |>
    mutate(
      GeneType = case_when(
        FlyBaseID %in% ssdo_ids ~ "SSDOs",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      ),
      GeneType = factor(GeneType, levels = ordered_gene_types)
    )
}


assign_collapsed_resolution_groups <- function(data) {
  data |>
    mutate(
      CollapsedGeneType = ifelse(FlyBaseID %in% ssdo_ids, "SSDOs", "Background"),
      CollapsedGeneType = factor(CollapsedGeneType, levels = c("SSDOs", "Background"))
    )
}


load_neal_candidate_ids <- function() {
  read.csv(
    input_path("gene_hits_snps_90.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(ensembl_gene_id)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct() |>
    pull(FlyBaseID)
}


load_ruzicka_candidate_ids <- function() {
  read.csv(
    input_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    distinct() |>
    pull(FlyBaseID)
}


assign_candidate_groups <- function(data, candidate_ids, group_var = "CandidateGroup") {
  data |>
    mutate(
      !!group_var := ifelse(
        FlyBaseID %in% candidate_ids,
        "Candidate SA genes",
        "Background genes"
      )
    )
}


assign_split_candidate_groups <- function(data, candidate_ids, group_var = "CandidateGroup") {
  disease_only_ids <- setdiff(disease_ids, ssdo_ids)
  ortholog_only_ids <- setdiff(ortholog_ids, union(ssdo_ids, disease_ids))

  data |>
    mutate(
      !!group_var := case_when(
        FlyBaseID %in% candidate_ids ~ "Candidate SA genes",
        FlyBaseID %in% disease_only_ids ~ "Disease orthologues",
        FlyBaseID %in% ortholog_only_ids ~ "Orthologues",
        TRUE ~ "Background"
      )
    )
}


run_candidate_resolution_enrichment <- function(universe_ids, focal_ids, feature_ids, study_name, feature_name, source_note) {
  analysis_df <- data.frame(FlyBaseID = unique(universe_ids), stringsAsFactors = FALSE) |>
    mutate(
      IsCandidate = FlyBaseID %in% focal_ids,
      IsFeature = FlyBaseID %in% feature_ids
    )

  contingency <- table(analysis_df$IsCandidate, analysis_df$IsFeature)
  fisher_result <- fisher.test(contingency)

  candidate_feature_count <- sum(analysis_df$IsCandidate & analysis_df$IsFeature)
  background_feature_count <- sum(!analysis_df$IsCandidate & analysis_df$IsFeature)
  candidate_count <- sum(analysis_df$IsCandidate)
  background_count <- sum(!analysis_df$IsCandidate)

  data.frame(
    Study = study_name,
    Feature = feature_name,
    Source = source_note,
    UniverseSize = nrow(analysis_df),
    CandidateCount = candidate_count,
    BackgroundCount = background_count,
    FeatureCount = sum(analysis_df$IsFeature),
    CandidateFeatureCount = candidate_feature_count,
    CandidateFeaturePercent = 100 * candidate_feature_count / candidate_count,
    BackgroundFeatureCount = background_feature_count,
    BackgroundFeaturePercent = 100 * background_feature_count / background_count,
    OddsRatio = unname(fisher_result$estimate),
    PValue = fisher_result$p.value,
    ConfLow = fisher_result$conf.int[1],
    ConfHigh = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )
}


run_candidate_enrichment_family <- function(universe_ids, focal_ids, study_name) {
  bind_rows(lapply(reported_resolution_features, function(feature_name) {
    source_note <- resolution_notes$Source[match(feature_name, resolution_notes$Feature)]

    run_candidate_resolution_enrichment(
      universe_ids = universe_ids,
      focal_ids = focal_ids,
      feature_ids = get_resolution_candidate_ids(feature_name),
      study_name = study_name,
      feature_name = feature_name,
      source_note = source_note
    )
  })) |>
    mutate(PAdjustBH = p.adjust(PValue, method = "BH")) |>
    arrange(PValue, Feature)
}


build_resolution_radar <- function(data, feature_order, feature_labels, group_var, group_order, group_colours, legend_ncol = 2) {
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("Install the 'fmsb' package to draw the SA-resolution radar chart.", call. = FALSE)
  }

  missing_features <- setdiff(feature_order, colnames(data))
  if (length(missing_features) > 0) {
    stop(
      sprintf("Radar plot features missing from sa_resolution_feature_membership.csv: %s", paste(missing_features, collapse = ", ")),
      call. = FALSE
    )
  }

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
      stringsAsFactors = FALSE
    )
  }))

  max_percent <- ceiling(max(radar_summary$Percent, na.rm = TRUE) / 10) * 10
  radial_step <- if (max_percent <= 50) 10 else 20
  radial_breaks <- seq(0, max_percent, by = radial_step)

  radar_matrix <- rbind(
    rep(max_percent, length(feature_order)),
    rep(0, length(feature_order)),
    do.call(
      rbind,
      lapply(group_order, function(group_name) {
        radar_summary |>
          filter(Group == group_name) |>
          arrange(match(Feature, feature_order)) |>
          pull(Percent)
      })
    )
  )

  radar_df <- as.data.frame(radar_matrix, stringsAsFactors = FALSE)
  colnames(radar_df) <- unname(feature_labels[feature_order])
  rownames(radar_df) <- c("max", "min", group_order)

  plot_radar <- function(output_path) {
    png(filename = output_path, width = 2700, height = 2400, res = 300, bg = "white")
    op <- par(family = "sans", mar = c(3.8, 2.8, 3.2, 2.8), xpd = NA)
    on.exit({
      par(op)
      dev.off()
    }, add = TRUE)

    fmsb::radarchart(
      radar_df,
      axistype = 1,
      pcol = group_colours[group_order],
      pfcol = grDevices::adjustcolor(group_colours[group_order], alpha.f = 0.20),
      plwd = 4,
      plty = 1,
      cglcol = "grey78",
      cglty = 1,
      axislabcol = "grey45",
      caxislabels = as.character(radial_breaks),
      cglwd = 0.8,
      calcex = 1.35,
      vlcex = 1.4
    )

    legend(
      "bottom",
      legend = group_order,
      horiz = FALSE,
      bty = "n",
      ncol = legend_ncol,
      lty = 1,
      lwd = 5,
      col = group_colours[group_order],
      cex = 1.25,
      pt.cex = 1.4,
      seg.len = 3,
      x.intersp = 0.8,
      y.intersp = 1.1,
      inset = c(0, -0.08)
    )
  }

  list(
    draw = plot_radar,
    summary = radar_summary |>
      select(Group, Feature, FeatureLabel, Percent)
  )
}


build_resolution_radar_permutation_tests <- function(data, feature_order, group_var, group_order, seed = 1L) {
  group_pairs <- combn(group_order, 2, simplify = FALSE)

  tests <- bind_rows(lapply(seq_along(feature_order), function(feature_index) {
    feature_name <- feature_order[[feature_index]]

    bind_rows(lapply(seq_along(group_pairs), function(pair_index) {
      pair <- group_pairs[[pair_index]]
      group_1 <- pair[[1]]
      group_2 <- pair[[2]]

      values_1 <- data[data[[group_var]] == group_1, feature_name, drop = TRUE]
      values_2 <- data[data[[group_var]] == group_2, feature_name, drop = TRUE]

      test_result <- two_sample_permutation_test(
        x = values_1,
        y = values_2,
        seed = seed + feature_index * 100 + pair_index
      )

      data.frame(
        Feature = feature_name,
        group_1 = group_1,
        group_2 = group_2,
        group_1_mean = mean(values_1, na.rm = TRUE),
        group_2_mean = mean(values_2, na.rm = TRUE),
        stringsAsFactors = FALSE
      ) |>
        bind_cols(test_result)
    }))
  }))

  tests |>
    group_by(Feature) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) |>
    ungroup()
}


whole_genome_universe <- resolution_membership |>
  pull(FlyBaseID) |>
  unique()

whole_genome_family <- run_enrichment_family(
  universe_ids = whole_genome_universe,
  universe_type_label = "Whole protein-coding fly genome"
)

resolution_plot_data <- assign_resolution_groups(resolution_membership)
collapsed_resolution_plot_data <- assign_collapsed_resolution_groups(resolution_membership)
neal_candidate_ids <- intersect(load_neal_candidate_ids(), resolution_membership$FlyBaseID)
ruzicka_candidate_ids <- intersect(load_ruzicka_candidate_ids(), resolution_membership$FlyBaseID)
neal_candidate_plot_data <- assign_candidate_groups(
  resolution_membership,
  neal_candidate_ids,
  group_var = "NealCandidateGroup"
)
ruzicka_candidate_plot_data <- assign_candidate_groups(
  resolution_membership,
  ruzicka_candidate_ids,
  group_var = "RuzickaCandidateGroup"
)
enrichment_results <- whole_genome_family$results


write.csv(
  enrichment_results,
  summary_path("SA_resolution_enrichment_tests.csv"),
  row.names = FALSE
)

radar_features <- c("HFB", "FB", "MB", "HMB", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "DR")
radar_feature_labels <- c(
  "HFB" = "HFB",
  "FB" = "FB",
  "MB" = "MB",
  "HMB" = "HMB",
  "Dup" = "Dup",
  "SSS" = "SSS",
  "ASE" = "ASE",
  "SS_ASE" = "SS-ASE",
  "TS_ASE" = "TS-ASE",
  "DR" = "DR"
)

radar_result <- build_resolution_radar(
  data = resolution_plot_data,
  feature_order = radar_features,
  feature_labels = radar_feature_labels,
  group_var = "GeneType",
  group_order = ordered_gene_types,
  group_colours = gene_type_palette,
  legend_ncol = 2
)
radar_permutation_tests <- build_resolution_radar_permutation_tests(
  data = resolution_plot_data,
  feature_order = radar_features,
  group_var = "GeneType",
  group_order = ordered_gene_types
)

write.csv(
  radar_result$summary,
  summary_path("SA_resolution_radar_prevalence.csv"),
  row.names = FALSE
)

write.csv(
  radar_permutation_tests,
  summary_path("SA_resolution_radar_permutation_tests.csv"),
  row.names = FALSE
)

radar_result$draw(figure_path("Fig4_sa_res_radar.png"))


# The current `SSS` column comes from `dimorphic.subset.list.txt`, treated here
# as the final consistently sex-specific splicing set from MaleLimitedEvo.
# Current documentation indicates that these 2036 genes passed cross-dataset
# consistency filtering for sex-specific splicing profiles.
