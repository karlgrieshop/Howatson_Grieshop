#!/usr/bin/env Rscript

# Build Table S2 as an editable Word document in a Table 1-style layout using
# overlap-aware candidate-context permutation tests for Fig. S3.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "Data")) &&
      dir.exists(file.path(current, "Sup_plotting_and_analyses")) &&
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
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
table_path <- function(...) file.path(project_root, "Tables", ...)
source(file.path(project_root, "Data_wrangling", "load_grieshop_candidate_ids.R"))


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(officer)
  library(flextable)
})


dir.create(summary_path(), recursive = TRUE, showWarnings = FALSE)
dir.create(table_path(), recursive = TRUE, showWarnings = FALSE)


ordered_groups <- c(
  "SSDOs",
  "Innocenti and Morrow candidates",
  "Ruzicka et al. candidates",
  "Grieshop et al. (2025) candidates"
)


normalise_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  trimws(x)
}


map_snp_metric_to_genes <- function(snp_values, gene_reference, source_label) {
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

  snps <- snp_values |>
    transmute(
      chromosome = normalise_chromosome(chromosome),
      start = as.integer(position),
      end = as.integer(position),
      Value = as.numeric(Value)
    ) |>
    filter(
      !is.na(chromosome),
      nzchar(chromosome),
      !is.na(start),
      !is.na(Value),
      is.finite(Value)
    )

  genes_dt <- as.data.table(genes)
  snps_dt <- as.data.table(snps)

  setkey(genes_dt, chromosome, start, end)
  setkey(snps_dt, chromosome, start, end)

  overlaps <- foverlaps(snps_dt, genes_dt, nomatch = 0L)

  as_tibble(overlaps) |>
    filter(!is.na(FlyBaseID)) |>
    group_by(FlyBaseID) |>
    summarise(
      Value = max(Value, na.rm = TRUE),
      Source = source_label,
      .groups = "drop"
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


load_innocenti_candidate_ids <- function(gene_reference) {
  read.csv(
    input_path("Innocenti_and_Morrow_ANTAGONISTIC.csv"),
    stringsAsFactors = FALSE
  ) |>
    transmute(Symbol = as.character(Symbol)) |>
    filter(!is.na(Symbol), nzchar(Symbol)) |>
    distinct() |>
    pull(Symbol) |>
    map_symbols_to_flybase(gene_reference = gene_reference) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    pull(FlyBaseID) |>
    unique()
}


load_ruzicka_candidate_ids <- function() {
  read.csv(
    input_path("Ruzicka_et_al_pbio.3000244.s008.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(FlyBaseID = as.character(`Flybase gene ID`)) |>
    filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID)) |>
    pull(FlyBaseID) |>
    unique()
}


build_ruzicka_magnitude_values <- function(gene_reference) {
  assoc <- read.delim(
    input_path("Assoc_Data.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    transmute(
      chromosome = as.character(Chromosome),
      position = as.integer(Basepair),
      Value = abs(as.numeric(Wald_Stat))
    ) |>
    filter(!is.na(Value), is.finite(Value))

  map_snp_metric_to_genes(
    snp_values = assoc[, c("chromosome", "position", "Value")],
    gene_reference = gene_reference,
    source_label = "Ruzicka et al. SA effect (max(abs(Wald_Stat)))"
  )
}


overlap_aware_two_sample_permutation_test <- function(
    value_df,
    group_1_ids,
    group_2_ids,
    n_perm = 10000,
    seed = 1L
) {
  value_df <- value_df |>
    transmute(FlyBaseID = as.character(FlyBaseID), Value = as.numeric(Value)) |>
    filter(!is.na(FlyBaseID), !is.na(Value), is.finite(Value)) |>
    distinct(FlyBaseID, .keep_all = TRUE)

  value_map <- setNames(value_df$Value, value_df$FlyBaseID)
  available_ids <- value_df$FlyBaseID

  group_1_ids <- intersect(unique(as.character(group_1_ids)), available_ids)
  group_2_ids <- intersect(unique(as.character(group_2_ids)), available_ids)
  shared_ids <- intersect(group_1_ids, group_2_ids)
  group_1_only_ids <- setdiff(group_1_ids, shared_ids)
  group_2_only_ids <- setdiff(group_2_ids, shared_ids)

  shared_values <- unname(value_map[shared_ids])
  group_1_only_values <- unname(value_map[group_1_only_ids])
  group_2_only_values <- unname(value_map[group_2_only_ids])

  observed_group_1 <- c(shared_values, group_1_only_values)
  observed_group_2 <- c(shared_values, group_2_only_values)

  if (length(observed_group_1) < 2 || length(observed_group_2) < 2) {
    return(data.frame(
      obs_diff = NA_real_,
      group_1_mean = mean(observed_group_1),
      group_2_mean = mean(observed_group_2),
      n_group_1 = length(observed_group_1),
      n_group_2 = length(observed_group_2),
      n_shared = length(shared_values),
      n_group_1_only = length(group_1_only_values),
      n_group_2_only = length(group_2_only_values),
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  obs_diff <- mean(observed_group_1) - mean(observed_group_2)

  if (length(group_1_only_values) == 0L && length(group_2_only_values) == 0L) {
    return(data.frame(
      obs_diff = obs_diff,
      group_1_mean = mean(observed_group_1),
      group_2_mean = mean(observed_group_2),
      n_group_1 = length(observed_group_1),
      n_group_2 = length(observed_group_2),
      n_shared = length(shared_values),
      n_group_1_only = 0L,
      n_group_2_only = 0L,
      p_value = 1,
      stringsAsFactors = FALSE
    ))
  }

  pooled_exclusive <- c(group_1_only_values, group_2_only_values)
  n_group_1_only <- length(group_1_only_values)

  set.seed(seed)
  perm_diffs <- replicate(n_perm, {
    shuffled <- sample(pooled_exclusive, size = length(pooled_exclusive), replace = FALSE)
    perm_group_1_only <- if (n_group_1_only > 0) shuffled[seq_len(n_group_1_only)] else numeric(0)
    perm_group_2_only <- if (n_group_1_only < length(shuffled)) shuffled[-seq_len(n_group_1_only)] else numeric(0)

    mean(c(shared_values, perm_group_1_only)) - mean(c(shared_values, perm_group_2_only))
  })

  data.frame(
    obs_diff = obs_diff,
    group_1_mean = mean(observed_group_1),
    group_2_mean = mean(observed_group_2),
    n_group_1 = length(observed_group_1),
    n_group_2 = length(observed_group_2),
    n_shared = length(shared_values),
    n_group_1_only = length(group_1_only_values),
    n_group_2_only = length(group_2_only_values),
    p_value = (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1),
    stringsAsFactors = FALSE
  )
}


pairwise_tests_with_overlap <- function(value_df, group_ids, group_order, seed = 1L) {
  group_pairs <- combn(group_order, 2, simplify = FALSE)

  tests <- bind_rows(lapply(seq_along(group_pairs), function(i) {
    group_pair <- group_pairs[[i]]
    group_1 <- group_pair[[1]]
    group_2 <- group_pair[[2]]

    overlap_aware_two_sample_permutation_test(
      value_df = value_df,
      group_1_ids = group_ids[[group_1]],
      group_2_ids = group_ids[[group_2]],
      seed = seed + i
    ) |>
      mutate(group_1 = group_1, group_2 = group_2, .before = 1)
  }))

  tests |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))
}


format_p_value <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  if (x > 0 && x < 0.0001) {
    return("<0.0001")
  }
  sprintf("%.4f", x)
}


group_labels <- c(
  "SSDOs" = "SSDOs",
  "Innocenti and Morrow candidates" = "Innocenti and Morrow candidates",
  "Grieshop et al. (2025) candidates" = "Grieshop et al. (2025) candidates",
  "Ruzicka et al. candidates" = "Ruzicka et al. candidates"
)


metric_comparisons <- list(
  "Ruzicka SA effect" = c(
    "SSDOs vs Innocenti and Morrow candidates",
    "SSDOs vs Grieshop et al. (2025) candidates",
    "SSDOs vs Ruzicka et al. candidates"
  )
)


gene_reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
protein_coding_gene_reference <- gene_reference |>
  filter(!is.na(FlyBaseID), grepl("^FBgn", FlyBaseID), gene_biotype == "protein_coding")

ssdo_ids <- read.csv(output_path("sexspecific_disease_ids.csv"), stringsAsFactors = FALSE) |>
  pull(FlyBaseID) |>
  unique()
innocenti_candidate_ids <- intersect(load_innocenti_candidate_ids(gene_reference), protein_coding_gene_reference$FlyBaseID)
ruzicka_candidate_ids <- intersect(load_ruzicka_candidate_ids(), protein_coding_gene_reference$FlyBaseID)
grieshop_candidate_ids <- intersect(
  load_grieshop_candidate_ids(input_path("DE.summary.tsv")),
  protein_coding_gene_reference$FlyBaseID
)

group_ids <- list(
  "SSDOs" = ssdo_ids,
  "Innocenti and Morrow candidates" = innocenti_candidate_ids,
  "Grieshop et al. (2025) candidates" = grieshop_candidate_ids,
  "Ruzicka et al. candidates" = ruzicka_candidate_ids
)

metric_inputs <- list(
  "Ruzicka SA effect" = build_ruzicka_magnitude_values(protein_coding_gene_reference)
)

metric_seeds <- c(
  "Ruzicka SA effect" = 300L
)

all_results <- bind_rows(lapply(names(metric_inputs), function(metric_name) {
  pairwise_tests_with_overlap(
    value_df = metric_inputs[[metric_name]],
    group_ids = group_ids,
    group_order = ordered_groups,
    seed = metric_seeds[[metric_name]]
  ) |>
    mutate(
      Metric = metric_name,
      Comparison = paste(group_labels[group_1], "vs", group_labels[group_2]),
      .before = 1
    )
}))

detailed_results <- bind_rows(lapply(names(metric_comparisons), function(metric_name) {
  all_results |>
    filter(Metric == metric_name, Comparison %in% metric_comparisons[[metric_name]]) |>
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))
}))

write.csv(
  detailed_results,
  summary_path("TableS2_candidate_context_permutation_tests.csv"),
  row.names = FALSE
)


display_comparison_order <- c(
  "SSDOs vs Innocenti and Morrow candidates",
  "SSDOs vs Grieshop et al. (2025) candidates",
  "SSDOs vs Ruzicka et al. candidates"
)

get_result_row <- function(metric_name, comparison_label) {
  detailed_results |>
    filter(Metric == metric_name, Comparison == comparison_label)
}


build_metric_cells <- function(metric_name, comparison_label) {
  if (!(comparison_label %in% metric_comparisons[[metric_name]])) {
    return(c("", "", ""))
  }

  metric_row <- get_result_row(metric_name, comparison_label)

  c(
    sprintf("%.4f", metric_row$obs_diff[[1]]),
    format_p_value(metric_row$p_value[[1]]),
    format_p_value(metric_row$p_adjust_bh[[1]])
  )
}

build_block_df <- function(metric_name, block_label) {
  bind_rows(lapply(display_comparison_order, function(comparison_label) {
    metric_cells <- build_metric_cells(metric_name, comparison_label)

    data.frame(
      Comparison = comparison_label,
      Observed_Difference = metric_cells[[1]],
      P_Value = metric_cells[[2]],
      BH_Adjusted_P = metric_cells[[3]],
      Block = block_label,
      stringsAsFactors = FALSE
    )
  }))
}


magnitude_df <- build_block_df(
  metric_name = "Ruzicka SA effect",
  block_label = "SA effect"
)

table_s2 <- bind_rows(
  data.frame(
    Comparison = "",
    Observed_Difference = "SA effect",
    P_Value = "",
    BH_Adjusted_P = "",
    Block = "block",
    stringsAsFactors = FALSE
  ),
  magnitude_df |>
    mutate(Block = "data")
)

ft <- flextable(table_s2 |> select(-Block)) |>
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
  width(j = "Comparison", width = 3.5) |>
  width(j = "Observed_Difference", width = 0.9) |>
  width(j = "P_Value", width = 0.8) |>
  width(j = "BH_Adjusted_P", width = 1.1) |>
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

block_rows <- which(table_s2$Block == "block")
for (row_index in block_rows) {
  ft <- merge_at(ft, i = row_index, j = 2:4, part = "body")
}
ft <- bold(ft, i = block_rows, j = 2, bold = TRUE, part = "body")
ft <- align(ft, i = block_rows, j = 2, align = "center", part = "body")

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = table_path("TableS2.docx"))
