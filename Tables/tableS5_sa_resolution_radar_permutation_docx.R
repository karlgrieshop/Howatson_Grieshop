# Build Table S5 as an editable Word document in a compact manuscript-ready
# layout using the tracked SA-resolution radar permutation summary.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
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
table_path <- function(...) file.path(project_root, "Tables", ...)


suppressPackageStartupMessages({
  library(dplyr)
  library(officer)
  library(flextable)
})


input_file <- summary_path("SA_resolution_radar_permutation_tests.csv")
prevalence_file <- summary_path("SA_resolution_radar_prevalence.csv")
feature_notes_file <- summary_path("sa_resolution_feature_notes.csv")
if (!file.exists(input_file)) {
  stop("Missing required permutation-test table: SA_resolution_radar_permutation_tests.csv", call. = FALSE)
}
if (!file.exists(prevalence_file)) {
  stop("Missing required prevalence table: SA_resolution_radar_prevalence.csv", call. = FALSE)
}
if (!file.exists(feature_notes_file)) {
  stop("Missing required feature-notes table: sa_resolution_feature_notes.csv", call. = FALSE)
}

dir.create(table_path(), recursive = TRUE, showWarnings = FALSE)

portrait_section <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69, orient = "portrait"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55),
  type = "continuous"
)


format_p_raw <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 1e-04, "<0.0001", sprintf("%.4f", x)))
}

format_p_adjusted <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.4f", x))
}

comparison_order <- c(
  "SSDOs vs Disease orthologues",
  "SSDOs vs Orthologues",
  "SSDOs vs Background",
  "Disease orthologues vs Orthologues",
  "Disease orthologues vs Background",
  "Orthologues vs Background"
)

feature_order <- c("HFB", "FB", "MB", "HMB", "HTS", "TS", "TG", "HTG", "Dup", "SSS", "ASE", "SS_ASE", "TS_ASE", "SSDR")
feature_labels <- c(
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


format_percent <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2f%%", as.numeric(x)))
}

format_count <- function(x) {
  prettyNum(as.integer(x), big.mark = ",", preserve.width = "none")
}

prevalence_lookup <- read.csv(prevalence_file, stringsAsFactors = FALSE) |>
  transmute(
    Feature = as.character(Feature),
    Group = as.character(Group),
    PercentLabel = format_percent(as.numeric(Percent))
  )

feature_count_lookup <- read.csv(feature_notes_file, stringsAsFactors = FALSE) |>
  transmute(
    Feature = as.character(Feature),
    FeatureCountLabel = paste0(unname(feature_labels[as.character(Feature)]), " (n = ", format_count(as.numeric(GeneCount)), ")")
  )

radar_perm <- read.csv(input_file, stringsAsFactors = FALSE) |>
  transmute(
    Feature = as.character(Feature),
    group_1 = as.character(group_1),
    group_2 = as.character(group_2),
    Observed_Difference = sprintf("%.4f", as.numeric(obs_diff)),
    P_Value = format_p_raw(as.numeric(p_value)),
    BH_Adjusted_P = format_p_adjusted(as.numeric(p_adjust_bh))
  ) |>
  left_join(
    prevalence_lookup |>
      rename(group_1_percent = PercentLabel, group_1 = Group),
    by = c("Feature", "group_1")
  ) |>
  left_join(
    prevalence_lookup |>
      rename(group_2_percent = PercentLabel, group_2 = Group),
    by = c("Feature", "group_2")
  ) |>
  mutate(
    Comparison = paste0(
      group_1, " (", group_1_percent, ") vs ",
      group_2, " (", group_2_percent, ")"
    )
  ) |>
  mutate(
    Feature = factor(Feature, levels = feature_order),
    ComparisonOrder = factor(paste(group_1, "vs", group_2), levels = comparison_order)
  ) |>
  arrange(Feature, ComparisonOrder) |>
  mutate(
    Feature = as.character(Feature),
    FeatureLabel = unname(feature_labels[Feature])
  )

table_s5_blocks <- lapply(feature_order[feature_order %in% unique(radar_perm$Feature)], function(feature_name) {
  feature_label <- feature_count_lookup$FeatureCountLabel[match(feature_name, feature_count_lookup$Feature)]
  feature_rows <- radar_perm[radar_perm$Feature == feature_name, , drop = FALSE]

  bind_rows(
    data.frame(
      Comparison = "",
      Observed_Difference = feature_label,
      P_Value = "",
      BH_Adjusted_P = "",
      Block = "block",
      stringsAsFactors = FALSE
    ),
    feature_rows |>
      transmute(
        Comparison,
        Observed_Difference,
        P_Value,
        BH_Adjusted_P,
        Block = "data"
      )
  )
})

table_s5 <- bind_rows(table_s5_blocks)


ft <- flextable(table_s5 |> select(-Block)) |>
  set_header_labels(
    Comparison = "Comparison",
    Observed_Difference = "Observed\ndifference",
    P_Value = "P-value",
    BH_Adjusted_P = "BH-adjusted P-value"
  ) |>
  align(align = "center", part = "all") |>
  align(j = "Comparison", align = "left", part = "body") |>
  valign(part = "all", valign = "center") |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 10, part = "all") |>
  line_spacing(part = "all", space = 1) |>
  padding(padding = 2, part = "all") |>
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) |>
  border_inner_h(part = "all", border = fp_border(color = "black", width = 0.75)) |>
  border_inner_v(part = "all", border = fp_border(color = "black", width = 0.75)) |>
  bold(part = "header") |>
  width(j = "Comparison", width = 4.1) |>
  width(j = "Observed_Difference", width = 0.8) |>
  width(j = "P_Value", width = 0.7) |>
  width(j = "BH_Adjusted_P", width = 1.0) |>
  set_table_properties(layout = "fixed") |>
  fit_to_width(max_width = 6.8)

ft <- compose(
  ft,
  i = 1,
  j = "P_Value",
  part = "header",
  value = as_paragraph(as_i("P"), "-value")
)

ft <- compose(
  ft,
  i = 1,
  j = "BH_Adjusted_P",
  part = "header",
  value = as_paragraph("BH-adjusted ", as_i("P"), "-value")
)

block_rows <- which(table_s5$Block == "block")
for (row_index in block_rows) {
  ft <- merge_at(ft, i = row_index, j = 2:4, part = "body")
}
ft <- bold(ft, i = block_rows, j = 2, bold = TRUE, part = "body")
ft <- align(ft, i = block_rows, j = 2, align = "center", part = "body")


doc <- read_docx() |>
  body_set_default_section(value = portrait_section) |>
  body_add_flextable(ft) |>
  body_end_block_section(block_section(portrait_section))


print(doc, target = table_path("TableS5.docx"))
