# Core manuscript analyses for SA-resolution mechanisms.
# This script is broader than the Fig4 radar itself, it assembles the
# whole-genome SA-resolution enrichment tests, the per-group radar prevalence
# values, the pairwise follow-up permutation tests and the omnibus multivariate
# permutation test that summarises whether the four gene classes differ in their
# overall mechanism profiles.


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


# Convenience wrapper: every downstream test works from the same binary
# gene-by-feature matrix, so this helper standardises how we recover the
# FlyBase IDs assigned TRUE for each mechanism.

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


# Fisher enrichment here asks whether SSDOs are unusually common among genes
# carrying a given mechanism, relative to the whole protein-coding universe.

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


# Apply the same whole-genome enrichment logic to every mechanism, then
# correct the resulting family of enrichment tests with BH.

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


# Generic two-sample permutation test used throughout the repo.
# We permute gene labels, not values, so the observed statistic is the
# difference in group means under the null of exchangeable group membership.

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


# Convert the binary mechanism table into the mutually exclusive manuscript
# comparison groups used everywhere else in the project.

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



# A collapsed SSDO-versus-rest grouping is retained for simpler binary contrasts.
assign_collapsed_resolution_groups <- function(data) {
  data |>
    mutate(
      CollapsedGeneType = ifelse(FlyBaseID %in% ssdo_ids, "SSDOs", "Background"),
      CollapsedGeneType = factor(CollapsedGeneType, levels = c("SSDOs", "Background"))
    )
}


# These helpers recover the external candidate SA gene sets used for
# mechanism-level enrichment against independent SA studies.
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


# Build the values that are actually plotted in Fig4:
# for each gene group and each mechanism, calculate the percentage of genes
# scored TRUE in the binary mechanism table.
# Convert the binary mechanism table into the actual percentages plotted in
# Fig4, then wrap the `fmsb` plotting call so those same values can also be
# exported as a tidy summary table.
build_resolution_radar <- function(data, feature_order, feature_labels, group_var, group_order, group_colours, legend_ncol = 2) {
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("Install the 'fmsb' package to draw the SA-resolution radar chart.", call. = FALSE)
  }

  plot_group_order <- rev(group_order)

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
      pcol = group_colours[plot_group_order],
      pfcol = grDevices::adjustcolor(group_colours[plot_group_order], alpha.f = 0.20),
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

# These pairwise tests follow the omnibus comparison and 
# are conducted on the same binary (0/1) matrix used for the radar plot. 
# Consequently, group means represent prevalence (proportion of TRUE values).
# Differences between groups can be interpreted directly as shifts in prevalence.

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
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))
}



# This custom omnibus statistic is a MANOVA-like between-group sum of squares on
# the binary mechanism matrix. It is computationally much lighter than
# distance-matrix approaches on the full 13,986-gene dataset.
compute_ssb <- function(x, groups) {
  grand_mean <- colMeans(x)
  group_levels <- unique(as.character(groups))

  sum(vapply(group_levels, function(group_name) {
    subset_x <- x[groups == group_name, , drop = FALSE]
    n_g <- nrow(subset_x)
    mean_g <- colMeans(subset_x)
    n_g * sum((mean_g - grand_mean)^2)
  }, numeric(1)))
}


run_multivariate_permutation_test <- function(x, groups, n_perm = 1999, seed = 1L) {
  groups <- factor(groups, levels = ordered_gene_types)
  x <- as.matrix(x)

  sst <- sum(scale(x, scale = FALSE)^2)
  ssb <- compute_ssb(x, groups)
  ssw <- sst - ssb
  df_between <- nlevels(groups) - 1
  df_within <- nrow(x) - nlevels(groups)
  obs_f <- (ssb / df_between) / (ssw / df_within)
  obs_r2 <- ssb / sst

  set.seed(seed)
  perm_f <- replicate(n_perm, {
    perm_groups <- sample(groups, length(groups), replace = FALSE)
    perm_ssb <- compute_ssb(x, perm_groups)
    perm_ssw <- sst - perm_ssb
    (perm_ssb / df_between) / (perm_ssw / df_within)
  })

  data.frame(
    N = nrow(x),
    P = ncol(x),
    DfBetween = df_between,
    DfWithin = df_within,
    SSB = ssb,
    SSW = ssw,
    R2 = obs_r2,
    PseudoF = obs_f,
    PValue = (sum(perm_f >= obs_f) + 1) / (n_perm + 1),
    Permutations = n_perm,
    stringsAsFactors = FALSE
  )
}


# ---- Main analysis workflow -----------------------------------------------
# Assemble the universes, run the enrichment tests, calculate the radar values,
# and then export the tables/figure used in the manuscript.

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

radar_features <- c("HFB", "FB", "MB", "HMB", "HTS", "TS", "TG", "HTG", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "SSDR")
radar_feature_labels <- c(
  "HFB" = "HFB",
  "FB" = "FB",
  "MB" = "MB",
  "HMB" = "HMB",
  "HTS" = "HTS",
  "TS" = "TS",
  "TG" = "TG",
  "HTG" = "HTG",
  "Dup" = "Dup",
  "SSS" = "SSS",
  "ASE" = "ASE",
  "SS_ASE" = "SS-ASE",
  "TS_ASE" = "TS-ASE",
  "SSDR" = "SSDR"
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
omnibus_multivariate_test <- run_multivariate_permutation_test(
  x = resolution_plot_data[, radar_features, drop = FALSE],
  groups = resolution_plot_data$GeneType,
  seed = 41L
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

write.csv(
  omnibus_multivariate_test,
  summary_path("SA_resolution_multivariate_tests.csv"),
  row.names = FALSE
)

radar_result$draw(figure_path("Fig4_sa_res_radar.png"))


# The current `SSS` column comes from `dimorphic.subset.list.txt`, treated here
# as the final 2036-gene consistently sex-specific splicing set from
# Grieshop et al. (2025), derived from Mishra et al. (2024a) reference
# profiles and filtered using three additional Mishra et al. (2024b)
# populations. HTS, TS, TG, and HTG are built upstream in
# `build_sa_resolution_inputs.R` as the non-middle tau quintiles from the Singh
# and Agrawal whole-body tissue-specificity estimates.
