# Build Table S4 as an editable Word document in a compact manuscript-ready
# layout using the tracked pi_rN, DoS, and Tajima's D permutation summaries.

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


required_files <- c(
  "rlnd Permutation Tests.csv",
  "DoS Permutation Tests.csv",
  "tajD Permutation Tests.csv"
)

missing_files <- required_files[!file.exists(summary_path(required_files))]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required permutation-test tables:\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
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


read_perm_table <- function(filename) {
  read.csv(summary_path(filename), stringsAsFactors = FALSE) |>
    transmute(
      Comparison = paste(as.character(group_1), "vs", as.character(group_2)),
      Observed_Difference = sprintf("%.4f", as.numeric(obs_diff)),
      P_Value = format_p_raw(as.numeric(p_value)),
      BH_Adjusted_P = format_p_adjusted(as.numeric(p_adjust_bh))
    ) |>
    mutate(Comparison = factor(Comparison, levels = comparison_order)) |>
    arrange(Comparison) |>
    mutate(Comparison = as.character(Comparison))
}


build_block_df <- function(filename, block_label) {
  read_perm_table(filename) |>
    mutate(Block = block_label) |>
    select(Comparison, Observed_Difference, P_Value, BH_Adjusted_P, Block)
}


rlnd_df <- build_block_df(
  filename = "rlnd Permutation Tests.csv",
  block_label = "pi_rN"
)

dos_df <- build_block_df(
  filename = "DoS Permutation Tests.csv",
  block_label = "DoS"
)

tajd_df <- build_block_df(
  filename = "tajD Permutation Tests.csv",
  block_label = "Tajima's D"
)

table_s4 <- bind_rows(
  data.frame(
    Comparison = "",
    Observed_Difference = "pi_rN",
    P_Value = "",
    BH_Adjusted_P = "",
    Block = "block",
    stringsAsFactors = FALSE
  ),
  rlnd_df |>
    mutate(Block = "data"),
  data.frame(
    Comparison = "",
    Observed_Difference = "DoS",
    P_Value = "",
    BH_Adjusted_P = "",
    Block = "block",
    stringsAsFactors = FALSE
  ),
  dos_df |>
    mutate(Block = "data"),
  data.frame(
    Comparison = "",
    Observed_Difference = "Tajima's D",
    P_Value = "",
    BH_Adjusted_P = "",
    Block = "block",
    stringsAsFactors = FALSE
  ),
  tajd_df |>
    mutate(Block = "data")
)


ft <- flextable(table_s4 |> select(-Block)) |>
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
  width(j = "Comparison", width = 2.45) |>
  width(j = "Observed_Difference", width = 1.0) |>
  width(j = "P_Value", width = 0.9) |>
  width(j = "BH_Adjusted_P", width = 1.2) |>
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

block_rows <- which(table_s4$Block == "block")
for (row_index in block_rows) {
  ft <- merge_at(ft, i = row_index, j = 2:4, part = "body")
}
ft <- bold(ft, i = block_rows, j = 2, bold = TRUE, part = "body")
ft <- align(ft, i = block_rows, j = 2, align = "center", part = "body")


doc <- read_docx() |>
  body_set_default_section(value = portrait_section) |>
  body_add_flextable(ft) |>
  body_end_block_section(block_section(portrait_section))


print(doc, target = table_path("TableS4.docx"))
