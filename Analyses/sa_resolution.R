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
output_path <- function(...) file.path(project_root, "Data", "outputs", ...)
figure_path <- function(...) file.path(project_root, "Figures", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)


suppressPackageStartupMessages({
  library(dplyr)
})


# Run `Data_wrangling/build_sa_resolution_inputs.R` first. That script assembles
# the genome-wide feature-membership table and records which requested
# SA-resolution features are ready versus still blocked by missing external
# sources.

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
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
reported_resolution_features <- setdiff(available_resolution_features, "DS")

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


bh_star_label <- function(adjusted_p) {
  if (is.na(adjusted_p)) {
    return("")
  }
  if (adjusted_p < 0.001) {
    return("***")
  }
  if (adjusted_p < 0.01) {
    return("**")
  }
  if (adjusted_p < 0.05) {
    return("*")
  }
  ""
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


build_resolution_radar <- function(data, feature_order, feature_labels) {
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

  radar_summary <- bind_rows(lapply(c(FALSE, TRUE), function(is_ssdo) {
    subset_df <- data |>
      filter(IsSSDO == is_ssdo)

    data.frame(
      Group = if (is_ssdo) "SSDOs" else "Non-SSDOs",
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
  group_colours <- c("Non-SSDOs" = "#7A7A7A", "SSDOs" = "#D55E00")

  radar_matrix <- rbind(
    rep(max_percent, length(feature_order)),
    rep(0, length(feature_order)),
    radar_summary |>
      filter(Group == "Non-SSDOs") |>
      arrange(match(Feature, feature_order)) |>
      pull(Percent),
    radar_summary |>
      filter(Group == "SSDOs") |>
      arrange(match(Feature, feature_order)) |>
      pull(Percent)
  )

  radar_df <- as.data.frame(radar_matrix, stringsAsFactors = FALSE)
  colnames(radar_df) <- unname(feature_labels[feature_order])
  rownames(radar_df) <- c("max", "min", "Non-SSDOs", "SSDOs")

  plot_radar <- function(output_path) {
    png(filename = output_path, width = 2400, height = 2400, res = 300, bg = "white")
    op <- par(family = "sans", mar = c(2.5, 2.5, 3.4, 2.5), xpd = NA)
    on.exit({
      par(op)
      dev.off()
    }, add = TRUE)

    fmsb::radarchart(
      radar_df,
      axistype = 1,
      pcol = group_colours[c("Non-SSDOs", "SSDOs")],
      pfcol = c(
        grDevices::adjustcolor(group_colours["Non-SSDOs"], alpha.f = 0.24),
        grDevices::adjustcolor(group_colours["SSDOs"], alpha.f = 0.24)
      ),
      plwd = 4,
      plty = 1,
      cglcol = "grey78",
      cglty = 1,
      axislabcol = "grey45",
      caxislabels = as.character(radial_breaks),
      cglwd = 0.8,
      vlcex = 1.15
    )

    legend(
      "top",
      legend = c("Non-SSDOs", "SSDOs"),
      horiz = TRUE,
      bty = "n",
      lty = 1,
      lwd = 4,
      col = group_colours[c("Non-SSDOs", "SSDOs")],
      cex = 1.1,
      inset = c(0, -0.02)
    )
  }

  list(
    draw = plot_radar,
    summary = radar_summary |>
      select(Group, Feature, FeatureLabel, Percent)
  )
}


whole_genome_universe <- resolution_membership |>
  pull(FlyBaseID) |>
  unique()

whole_genome_family <- run_enrichment_family(
  universe_ids = whole_genome_universe,
  universe_type_label = "Whole protein-coding fly genome"
)

enrichment_results <- whole_genome_family$results
whole_genome_results <- enrichment_results

overlap_ids <- whole_genome_family$overlaps


write.csv(
  enrichment_results,
  summary_path("SA_resolution_enrichment_tests.csv"),
  row.names = FALSE
)

write.csv(
  overlap_ids,
  summary_path("SA_resolution_enrichment_overlap_ids.csv"),
  row.names = FALSE
)


radar_features <- c("SBGE", "Dup", "SSS", "Cis", "SS_Cis", "DR")
radar_feature_stars <- setNames(
  vapply(
    radar_features,
    function(feature_name) {
      bh_star_label(whole_genome_results$PAdjustBH[match(feature_name, whole_genome_results$Study)])
    },
    character(1)
  ),
  radar_features
)
radar_feature_labels <- c(
  "SBGE" = paste0("SBGE", radar_feature_stars["SBGE"]),
  "Dup" = paste0("Dup", radar_feature_stars["Dup"]),
  "SSS" = paste0("SSS", radar_feature_stars["SSS"]),
  "Cis" = paste0("Cis", radar_feature_stars["Cis"]),
  "SS_Cis" = paste0("SS-Cis", radar_feature_stars["SS_Cis"]),
  "DR" = paste0("DR", radar_feature_stars["DR"])
)

radar_result <- build_resolution_radar(
  resolution_membership,
  feature_order = radar_features,
  feature_labels = radar_feature_labels
)

write.csv(
  radar_result$summary,
  summary_path("SA_resolution_radar_prevalence.csv"),
  row.names = FALSE
)

radar_result$draw(figure_path("sa_resolution_radar.png"))


# The current `SSS` column comes from `dimorphic.subset.list.txt`, which is the
# best available on-disk proxy for the requested broader SSS set but should be
# checked later against Michelle Liu's original source if the expected count
# must be ~3704 rather than 2036.
