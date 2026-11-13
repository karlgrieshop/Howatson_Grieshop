#!/usr/bin/env Rscript
# Exercise the active SNP-to-gene loaders without running their plotting workflows.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

genes <- data.frame(
  FlyBaseID = sprintf("FBgn%07d", 1:5), chromosome = "2L",
  start = seq(10L, 50L, 10L), end = seq(19L, 59L, 10L),
  gene_biotype = "protein_coding"
)
fixture <- data.frame(
  Chromosome = "chr2L", Basepair = c(10L, 11L, 20L, 21L, 30L, 31L, 32L, 40L, 41L, 999L),
  Wald_Stat = c(-4, 1, -3, -2, 0, NA, Inf, NA, -Inf, 100)
)
expected <- setNames(c(4, 3, 0), genes$FlyBaseID[1:3])
scripts <- c(
  "Plotting_and_analyses/sa_magnitude_plots.R",
  "Plotting_and_analyses/figure1_sa_magnitude.R",
  "Plotting_and_analyses/figureS4_figure1_with_brigos_5class.R",
  "Sup_plotting_and_analyses/figureS2_sa_mag.R",
  "Sup_plotting_and_analyses/figureS3_context.R",
  "Tables/tableS2_candidate_context_docx.R",
  "Sup_plotting_and_analyses/figureS7_resolution_stages.R",
  "Sup_plotting_and_analyses/figureS9_conservation_adjusted_sa.R"
)

for (script in scripts) {
  expressions <- parse(script)
  for (direction in c(1, -1)) {
    env <- new.env(parent = globalenv())
    env$input_path <- function(...) paste(..., collapse = "/")
    env$read.delim <- function(...) transform(fixture, Wald_Stat = direction * Wald_Stat)
    env$fread <- env$read.delim
    env$genes <- select(genes, FlyBaseID, chromosome, start, end)
    env$assign_gene_types <- identity
    function_names <- c(
      "normalise_chromosome", "map_snp_metric_to_genes",
      "build_ruzicka_values", "build_ruzicka_magnitude_values"
    )
    for (expression in expressions) {
      if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
          is.symbol(expression[[2]]) &&
          as.character(expression[[2]]) %in% function_names) {
        eval(expression, env)
      }
    }
    if (exists("build_ruzicka_values", envir = env, inherits = FALSE)) {
      result <- env$build_ruzicka_values(genes)
    } else if (exists("build_ruzicka_magnitude_values", envir = env, inherits = FALSE)) {
      result <- env$build_ruzicka_magnitude_values(genes)
    } else {
      # S7/S9 perform their mapping at top level; evaluate only that mapping.
      for (expression in expressions) {
        if (!is.call(expression)) next
        is_assignment <- identical(expression[[1]], as.name("<-"))
        if ((is_assignment && is.symbol(expression[[2]]) &&
             as.character(expression[[2]]) %in%
               c("snps", "genes_dt", "snps_dt", "gene_values", "sa_values")) ||
            identical(expression[[1]], as.name("setkey"))) {
          eval(expression, env)
        }
      }
      result <- if (exists("sa_values", envir = env, inherits = FALSE)) {
        rename(env$sa_values, Value = SA)
      } else {
        env$gene_values
      }
    }
    observed <- setNames(result$Value, result$FlyBaseID)
    if (!identical(sort(names(observed)), sort(names(expected))) ||
        !isTRUE(all.equal(unname(observed[names(expected)]), unname(expected)))) {
      stop("Incorrect maximum absolute Wald statistic in ", script,
           " with allele direction ", direction, call. = FALSE)
    }
  }
}
message("All eight SA loaders retain max(|z|), exclude nonfinite/unmapped SNPs, and are allele-flip invariant.")
