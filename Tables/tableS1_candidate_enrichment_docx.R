# Build Table S1 as an editable Word document from the candidate-SA Fisher
# enrichment results table.

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


input_file <- summary_path("TableS1_SA_Candidate_Enrichment_Tests.csv")
if (!file.exists(input_file)) {
  stop("Missing required input table: TableS1_SA_Candidate_Enrichment_Tests.csv", call. = FALSE)
}

dir.create(table_path(), recursive = TRUE, showWarnings = FALSE)

portrait_section <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69, orient = "portrait"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55),
  type = "continuous"
)


format_p <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 1e-04, "<0.0001", sprintf("%.3f", x)))
}

format_3dp <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.3f", x))
}


study_labels <- c(
  "Innocenti and Morrow" = "I&M",
  "Ruzicka et al." = "Ruz.",
  "Brigos-Baril et al. (2026)" = "Brigos et al."
)

focal_labels <- c(
  "SSDOs" = "SSDOs",
  "Disease orthologues" = "Dis. orth.",
  "Orthologues" = "Orth."
)


table_s1 <- read.csv(input_file, stringsAsFactors = FALSE) |>
  transmute(
    Focal_Group = recode(FocalGroup, !!!focal_labels),
    Study = recode(Study, !!!study_labels),
    Total = as.integer(UniverseSize),
    Focal_Count = as.integer(FocalCount),
    Candidate_Count = as.integer(CandidateCount),
    Overlap_count = as.integer(OverlapCount),
    Odds_Ratio = format_3dp(as.numeric(OddsRatio)),
    P_Value = format_p(as.numeric(PValue)),
    P_BH = format_p(as.numeric(PAdjustedBH)),
    LCI = format_3dp(as.numeric(ConfLow)),
    UCI = format_3dp(as.numeric(ConfHigh))
  ) |>
  mutate(across(everything(), as.character))

ft <- flextable(table_s1) |>
  theme_booktabs() |>
  bold(part = "header") |>
  set_header_labels(
    Focal_Group = "Focal\nGroup",
    Study = "Study",
    Total = "Total",
    Focal_Count = "Focal\nCount",
    Candidate_Count = "Cand.\nCount",
    Overlap_count = "Overlap\ncount",
    Odds_Ratio = "Odds\nRatio",
    P_Value = "",
    P_BH = "",
    LCI = "LCI",
    UCI = "UCI"
  ) |>
  align(align = "center", part = "all") |>
  align(j = c("Focal_Group", "Study"), align = "left", part = "body") |>
  valign(part = "all", valign = "center") |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 10, part = "all") |>
  line_spacing(part = "all", space = 1) |>
  padding(padding = 2, part = "all") |>
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) |>
  border_inner_h(part = "all", border = fp_border(color = "black", width = 0.75)) |>
  border_inner_v(part = "all", border = fp_border(color = "black", width = 0)) |>
  width(j = "Focal_Group", width = 0.85) |>
  width(j = "Study", width = 0.55) |>
  width(j = "Total", width = 0.50) |>
  width(j = "Focal_Count", width = 0.62) |>
  width(j = "Candidate_Count", width = 0.62) |>
  width(j = "Overlap_count", width = 0.62) |>
  width(j = "Odds_Ratio", width = 0.62) |>
  width(j = "P_Value", width = 0.62) |>
  width(j = "P_BH", width = 0.62) |>
  width(j = c("LCI", "UCI"), width = 0.56) |>
  set_table_properties(layout = "fixed")

ft <- compose(
  ft,
  i = 1,
  part = "header",
  j = "P_Value",
  value = as_paragraph(as_i("P"), "-value")
)
ft <- compose(
  ft,
  i = 1,
  part = "header",
  j = "P_BH",
  value = as_paragraph(as_i("P"), " (BH)")
)
ft <- add_footer_lines(
  ft,
  values = "Raw two-sided Fisher P-values and BH-adjusted P-values across all nine focal-group by candidate-study tests."
)

doc <- read_docx() |>
  body_set_default_section(value = portrait_section) |>
  body_add_flextable(ft) |>
  body_end_block_section(block_section(portrait_section))


print(doc, target = table_path("tableS1.docx"))
