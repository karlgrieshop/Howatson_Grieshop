find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Data_wrangling")) &&
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
mishra_ase_path <- function(...) file.path(project_root, "..", "Mishra_et_al_Drosophila_ASE", ...)


suppressPackageStartupMessages({
  library(dplyr)
})


input_file <- mishra_ase_path(
  "Data_for_R_scripts",
  "DESeq2.outputs",
  "DifferentialGeneExpression.whole.bodies.tsv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Missing Mishra whole-body DE file:\n",
      input_file
    ),
    call. = FALSE
  )
}


dir.create(output_path(), recursive = TRUE, showWarnings = FALSE)


x1 <- 1
x2 <- 5
y0 <- 1

mishra <- read.delim(input_file, stringsAsFactors = FALSE) |>
  transmute(
    FlyBaseID = as.character(FlyBaseID),
    log2FoldChange = as.numeric(log2FoldChange),
    lfcSE = as.numeric(lfcSE)
  ) |>
  filter(
    !is.na(FlyBaseID),
    grepl("^FBgn", FlyBaseID),
    !is.na(log2FoldChange),
    is.finite(log2FoldChange),
    !is.na(lfcSE),
    is.finite(lfcSE)
  ) |>
  distinct(FlyBaseID, .keep_all = TRUE)


high_female <- mishra |>
  filter(log2FoldChange < -x2, (log2FoldChange + lfcSE) < -(x2 - y0)) |>
  mutate(Sex.Bias.Category = "Highly female-biased")

female <- mishra |>
  filter(
    log2FoldChange < -x1,
    log2FoldChange > -x2,
    (log2FoldChange + lfcSE) < 0,
    (log2FoldChange - lfcSE) > -(x2 + y0)
  ) |>
  mutate(Sex.Bias.Category = "Moderately female-biased")

unbiased <- mishra |>
  filter(
    log2FoldChange < x1,
    log2FoldChange > -x1,
    (log2FoldChange - lfcSE) > -(x1 + y0),
    (log2FoldChange + lfcSE) < (x1 + y0)
  ) |>
  mutate(Sex.Bias.Category = "Unbiased")

male <- mishra |>
  filter(
    log2FoldChange > x1,
    log2FoldChange < x2,
    (log2FoldChange - lfcSE) > 0,
    (log2FoldChange + lfcSE) < (x2 + y0)
  ) |>
  mutate(Sex.Bias.Category = "Moderately male-biased")

high_male <- mishra |>
  filter(log2FoldChange > x2, (log2FoldChange - lfcSE) > (x2 - y0)) |>
  mutate(Sex.Bias.Category = "Highly male-biased")


sbge_bins <- bind_rows(
  high_female,
  female,
  unbiased,
  male,
  high_male
) |>
  distinct(FlyBaseID, .keep_all = TRUE) |>
  arrange(
    factor(
      Sex.Bias.Category,
      levels = c(
        "Highly female-biased",
        "Moderately female-biased",
        "Unbiased",
        "Moderately male-biased",
        "Highly male-biased"
      )
    ),
    FlyBaseID
  )


write.csv(
  sbge_bins,
  output_path("mishra_whole_body_sbge_bins.csv"),
  row.names = FALSE
)
