#!/usr/bin/env Rscript
# Build the editable Table S6 from the retained FigS7 permutation results.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "Data")) &&
        file.exists(file.path(current, "README.md"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Run this script from inside the repository.", call. = FALSE)
    }
    current <- parent
  }
}

project_root <- find_project_root()
input_file <- file.path(project_root, "Data", "Summary_tables",
                        "FigS7_pairwise_permutation_tests.csv")
output_file <- file.path(project_root, "Tables", "TableS6.docx")

suppressPackageStartupMessages({
  library(dplyr)
  library(officer)
  library(flextable)
})

results <- read.csv(input_file, stringsAsFactors = FALSE)
metrics <- c("identity_percent", "fly_coverage_percent", "human_coverage_percent")
metric_labels <- c("Whole-alignment identity", "Fly protein coverage", "Human protein coverage")
group_names <- c("SSDOs", "Disease orthologues excluding SSDOs",
                 "Orthologues excluding disease genes")
group_labels <- c("SSDOs", "Disease orthologues", "Orthologues")
expected_pairs <- data.frame(
  Group1 = group_names[c(1, 1, 2)],
  Group2 = group_names[c(2, 3, 3)],
  ContrastOrder = 1:3
)
required <- c("AnalysisSet", "Metric", "Group1", "Group2", "ContrastOrder",
              "NGroup1", "NGroup2", "MeanGroup1", "MeanGroup2",
              "MeanDifference_pp", "PValue", "PAdjustedBH")
if (!all(required %in% names(results)) || nrow(results) != 9L ||
    any(results$AnalysisSet != "All mapped orthologues") ||
    !setequal(results$Metric, metrics) || anyNA(results[required])) {
  stop("Expected nine complete tests for the three measured orthologue categories.", call. = FALSE)
}

for (metric in metrics) {
  tests <- results[results$Metric == metric, ] |>
    arrange(ContrastOrder)
  if (nrow(tests) != 3L ||
      !identical(tests[c("Group1", "Group2", "ContrastOrder")], expected_pairs) ||
      any(!is.finite(tests$MeanDifference_pp)) ||
      any(tests$PValue < 0 | tests$PValue > 1) ||
      !isTRUE(all.equal(tests$MeanDifference_pp, tests$MeanGroup1 - tests$MeanGroup2)) ||
      !isTRUE(all.equal(tests$PAdjustedBH, p.adjust(tests$PValue, "BH")))) {
    stop("Invalid comparisons, differences, or BH family for ", metric, ".", call. = FALSE)
  }
}
counts <- bind_rows(
  results |> transmute(Group = Group1, N = NGroup1),
  results |> transmute(Group = Group2, N = NGroup2)
) |>
  distinct()
if (nrow(counts) != 3L || !setequal(counts$Group, group_names) ||
    any(!is.finite(counts$N) | counts$N < 1)) {
  stop("Category sample sizes are inconsistent across comparisons.", call. = FALSE)
}
format_p_raw <- function(x) ifelse(x < 1e-04, "<0.0001", sprintf("%.4f", x))
blocks <- lapply(seq_along(metrics), function(i) {
  tests <- results |>
    filter(Metric == metrics[i]) |>
    arrange(ContrastOrder)
  bind_rows(
    data.frame(Comparison = "", Observed_Difference = metric_labels[i],
               P_Value = "", BH_Adjusted_P = "", Block = "block", Significant = FALSE),
    tests |>
      transmute(
        Comparison = paste(group_labels[match(Group1, group_names)], "vs",
                           group_labels[match(Group2, group_names)]),
        Observed_Difference = sprintf("%.4f", MeanDifference_pp),
        P_Value = format_p_raw(PValue),
        BH_Adjusted_P = sprintf("%.4f", PAdjustedBH),
        Block = "data",
        Significant = PAdjustedBH < 0.05
      )
  )
})
table_s6 <- bind_rows(blocks)

portrait_section <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69, orient = "portrait"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55),
  type = "continuous"
)
ft <- flextable(table_s6 |> select(-Block, -Significant)) |>
  set_header_labels(
    Comparison = "Comparison", Observed_Difference = "Observed\ndifference",
    P_Value = "P-value", BH_Adjusted_P = "BH-adjusted P-value"
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
  width(j = "Comparison", width = 2.75) |>
  width(j = "Observed_Difference", width = 1.0) |>
  width(j = "P_Value", width = 0.9) |>
  width(j = "BH_Adjusted_P", width = 1.2) |>
  set_table_properties(layout = "fixed") |>
  fit_to_width(max_width = 6.8)
ft <- compose(ft, i = 1, j = "P_Value", part = "header",
              value = as_paragraph(as_i("P"), "-value"))
ft <- compose(ft, i = 1, j = "BH_Adjusted_P", part = "header",
              value = as_paragraph("BH-adjusted ", as_i("P"), "-value"))
block_rows <- which(table_s6$Block == "block")
for (row in block_rows) ft <- merge_at(ft, i = row, j = 2:4, part = "body")
ft <- bold(ft, i = block_rows, j = 2, bold = TRUE, part = "body")
ft <- align(ft, i = block_rows, j = 2, align = "center", part = "body")
ft <- bold(ft, i = which(table_s6$Significant), j = "BH_Adjusted_P",
           bold = TRUE, part = "body")

doc <- read_docx() |>
  body_set_default_section(value = portrait_section) |>
  body_add_flextable(ft) |>
  body_end_block_section(block_section(portrait_section))
print(doc, target = output_file)
