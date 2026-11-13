#!/usr/bin/env Rscript
# Maximum absolute Wald z-statistic before/after conservation adjustment on the same genes.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(project_root, "Data", "outputs", "dmel_gene_reference.csv"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || arguments[1] != "--alignment-dir") {
  stop(paste(
    "Usage: Rscript Sup_plotting_and_analyses/figureS8_conservation_adjusted_sa.R",
    "--alignment-dir fly_human_conservation/alignment_<run-id>"
  ), call. = FALSE)
}
alignment_dir <- normalizePath(arguments[2], winslash = "/", mustWork = TRUE)
source(file.path(project_root, "Sup_plotting_and_analyses", "conservation_sa_regression.R"))
output_path <- function(...) file.path(project_root, "Data", "outputs", ...)
summary_path <- function(...) file.path(project_root, "Data", "Summary_tables", ...)
input_path <- function(...) file.path(project_root, "Data", "inputs", ...)
category_order <- c(
  "SSDOs", "Disease orthologues excluding SSDOs", "Orthologues excluding disease genes"
)
palette <- setNames(c("#D55E00", "#0072B2", "#009E73"), category_order)
covariates <- c("identity_percent", "fly_coverage_percent", "human_coverage_percent")
n_resamples <- 10000L
n_boot <- 5000L

reference <- read.csv(output_path("dmel_gene_reference.csv"), stringsAsFactors = FALSE)
ssdo <- read.csv(output_path("sexspecific_disease_ids.csv"))$FlyBaseID
disease <- read.csv(output_path("all_disease_ortholog_ids.csv"))$FlyBaseID
orthologues <- read.csv(output_path("all_ortholog_ids.csv"))$FlyBaseID
if (!all(ssdo %in% disease) || !all(disease %in% orthologues)) {
  stop("Canonical gene-category nesting is inconsistent.")
}
genes <- reference |>
  filter(gene_biotype == "protein_coding") |>
  transmute(
    FlyBaseID, chromosome = trimws(gsub("^chr", "", chromosome, ignore.case = TRUE)),
    start = as.integer(start), end = as.integer(end)
  )
if (anyDuplicated(genes$FlyBaseID)) stop("Duplicated protein-coding reference genes.")
if (anyNA(genes) || any(!nzchar(genes$chromosome)) || any(genes$start > genes$end)) {
  stop("Invalid protein-coding gene coordinates.")
}
message("Mapping Ruzicka max(abs(Wald_Stat)) SA effects to genes.")
snps <- fread(input_path("Assoc_Data.txt"), select = c("Chromosome", "Basepair", "Wald_Stat")) |>
  as_tibble() |>
  transmute(
    chromosome = trimws(gsub("^chr", "", as.character(Chromosome), ignore.case = TRUE)),
    start = as.integer(Basepair), end = as.integer(Basepair), SA = abs(as.numeric(Wald_Stat))
  ) |>
  filter(!is.na(chromosome), nzchar(chromosome), !is.na(start), is.finite(SA))
genes_dt <- as.data.table(genes)
snps_dt <- as.data.table(snps)
setkey(genes_dt, chromosome, start, end)
setkey(snps_dt, chromosome, start, end)
sa_values <- as_tibble(foverlaps(snps_dt, genes_dt, nomatch = 0L)) |>
  group_by(FlyBaseID) |>
  summarise(SA = max(SA), NSNPs = n(), .groups = "drop")
rm(snps, snps_dt, genes_dt)

conservation_file <- file.path(alignment_dir, "tables", "gene_conservation.csv")
conservation <- read.csv(conservation_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "FlyBaseID", "category", "human_gene_id", covariates
)
if (!all(required %in% names(conservation)) || anyDuplicated(conservation$FlyBaseID) ||
    any(!is.finite(as.matrix(conservation[covariates]))) ||
    any(as.matrix(conservation[covariates]) < 0 | as.matrix(conservation[covariates]) > 100) ||
    anyNA(conservation$human_gene_id) || any(!nzchar(conservation$human_gene_id))) {
  stop("Conservation inputs must be unique, finite percentages with human-partner IDs.")
}
audit <- genes |>
  select(FlyBaseID) |>
  mutate(Category = case_when(
    FlyBaseID %in% ssdo ~ category_order[1],
    FlyBaseID %in% disease ~ category_order[2],
    FlyBaseID %in% orthologues ~ category_order[3],
    TRUE ~ "Background"
  )) |>
  left_join(sa_values, by = "FlyBaseID") |>
  left_join(conservation |> select(all_of(required)), by = "FlyBaseID") |>
  mutate(
    HasSA = is.finite(SA),
    HasConservation = FlyBaseID %in% conservation$FlyBaseID,
    IncludedMain = HasSA & HasConservation & Category != "Background",
    ExclusionReason = case_when(
      Category == "Background" ~ "No defined fly-human orthologue conservation",
      !HasSA & !HasConservation ~ "No mapped SA statistic and no validated conservation",
      !HasSA ~ "No mapped SA statistic",
      !HasConservation ~ "No validated conservation",
      TRUE ~ "Included"
    )
  )
if (any(audit$Category[audit$HasConservation] != audit$category[audit$HasConservation]) ||
    !all(conservation$FlyBaseID %in% audit$FlyBaseID)) {
  stop("Conservation categories differ from the current canonical protein-coding universe.")
}
values <- audit |> filter(IncludedMain)
if (!setequal(unique(values$Category), category_order)) stop("A category has no measured SA values.")
values$Category <- factor(values$Category, levels = category_order)
coverage <- audit |>
  group_by(Category) |>
  summarise(
    ProteinCodingEligible = n(), WithSA = sum(HasSA),
    WithConservation = sum(HasConservation), IncludedMain = sum(IncludedMain),
    MissingSAWithConservation = sum(HasConservation & !HasSA),
    WithSAWithoutConservation = sum(HasSA & !HasConservation), .groups = "drop"
  )
write.csv(audit, summary_path("FigS8_gene_inclusion_audit.csv"), row.names = FALSE)
write.csv(coverage, summary_path("FigS8_category_coverage.csv"), row.names = FALSE)

model_design <- function(data, adjusted) {
  groups <- droplevels(data$Category)
  data$Category <- groups
  centers <- colMeans(data[covariates])
  scales <- vapply(data[covariates], sd, numeric(1))
  if (any(scales <= 0)) stop("A conservation covariate has zero variance.")
  for (metric in covariates) data[[metric]] <- (data[[metric]] - centers[metric]) / scales[metric]
  nuisance <- if (adjusted) reformulate(covariates) else ~1
  full_formula <- if (adjusted) reformulate(c(covariates, "Category")) else ~Category
  x <- model.matrix(full_formula, data)
  reduced <- model.matrix(nuisance, data)
  predictions <- data.frame(Category = factor(levels(groups), levels = levels(groups)))
  for (metric in covariates) predictions[[metric]] <- 0
  prediction <- model.matrix(full_formula, predictions)
  hypothesis <- matrix(0, nlevels(groups) - 1L, ncol(x))
  hypothesis[, (ncol(reduced) + 1L):ncol(x)] <- diag(nlevels(groups) - 1L)
  list(
    x = x, reduced = reduced, prediction = prediction, hypothesis = hypothesis,
    centers = centers, scales = scales, data = data
  )
}

run_set <- function(data, analysis_set, seed_offset) {
  tests <- list()
  summaries <- list()
  coefficients <- list()
  diagnostics <- list()
  fitted_values <- list()
  counter <- 0L
  for (adjusted in c(FALSE, TRUE)) {
    model_label <- if (adjusted) "Conservation adjusted" else "Unadjusted (same genes)"
    for (subset in c("All three orthologue categories", "Disease-only vs other orthologues")) {
      counter <- counter + 1L
      seed <- seed_offset + counter
      subset_data <- if (subset == "All three orthologue categories") {
        data
      } else {
        data |> filter(Category != "SSDOs")
      }
      design <- model_design(subset_data, adjusted)
      message(analysis_set, " / ", model_label, " / ", subset,
              ": ", n_resamples, " restricted wild-cluster residual resamples.")
      result <- wild_cluster_test(
        design$x, design$reduced, subset_data$SA, subset_data$human_gene_id,
        design$hypothesis, n_resamples, seed
      )
      prediction <- if (subset == "All three orthologue categories") {
        design$prediction
      } else {
        design$prediction[1, , drop = FALSE] - design$prediction[2, , drop = FALSE]
      }
      boot <- cluster_bootstrap_predictions(
        design$x, subset_data$SA, subset_data$human_gene_id,
        prediction, n_boot, seed + 100L
      )
      counts <- table(subset_data$Category)
      is_pair <- subset == "Disease-only vs other orthologues"
      r2 <- 1 - sum(result$fit$residual^2) /
        sum((subset_data$SA - mean(subset_data$SA))^2)
      tests[[counter]] <- data.frame(
        AnalysisSet = analysis_set, Model = model_label, Hypothesis = subset,
        NGenes = nrow(subset_data), NHumanClusters = result$fit$g,
        NSSDOs = unname(counts[1]), NDiseaseOnly = unname(counts[2]),
        NOtherOrthologues = unname(counts[3]),
        EstimateDiseaseMinusOther = if (is_pair) boot$estimate else NA_real_,
        DifferenceCI95Low = if (is_pair) boot$low else NA_real_,
        DifferenceCI95High = if (is_pair) boot$high else NA_real_,
        WaldF = result$statistic, HypothesisDF = nrow(design$hypothesis),
        PValue = result$p_value, Resamples = n_resamples, Seed = seed,
        BootstrapReplicates = n_boot, BootstrapSeed = seed + 100L,
        RSquared = r2,
        Covariates = if (adjusted) paste(covariates, collapse = ";") else "None",
        Method = paste(
          "Restricted wild-cluster bootstrap, Rademacher signs on null residuals;",
          "CR1-studentized Wald statistic; clusters = human_gene_id; plus-one P"
        ),
        stringsAsFactors = FALSE
      )
      coefficients[[counter]] <- data.frame(
        AnalysisSet = analysis_set, Model = model_label, Hypothesis = subset,
        Term = colnames(design$x), Estimate = result$fit$beta,
        ClusterRobustSE = sqrt(diag(result$fit$covariance)),
        NGenes = nrow(subset_data), NHumanClusters = result$fit$g
      )
      diagnostics[[counter]] <- data.frame(
        AnalysisSet = analysis_set, Model = model_label, Hypothesis = subset,
        Covariate = covariates, CenterPercent = unname(design$centers),
        ScalePercent = unname(design$scales), DesignConditionNumber = kappa(design$x),
        RSquared = r2,
        ResidualSD_SSDOs = sd(result$fit$residual[subset_data$Category == category_order[1]]),
        ResidualSD_DiseaseOnly = sd(result$fit$residual[subset_data$Category == category_order[2]]),
        ResidualSD_Other = sd(result$fit$residual[subset_data$Category == category_order[3]])
      )
      if (!is_pair) {
        null_fit <- fit_conservation_model(
          design$reduced, subset_data$SA, subset_data$human_gene_id
        )
        fitted_values[[model_label]] <- data.frame(
          AnalysisSet = analysis_set, Model = model_label,
          FlyBaseID = subset_data$FlyBaseID, Category = subset_data$Category,
          HumanGeneID = subset_data$human_gene_id, SA = subset_data$SA,
          FittedFullModel = result$fit$fitted, ResidualFullModel = result$fit$residual,
          FittedWithoutCategory = null_fit$fitted,
          ResidualWithoutCategory = null_fit$residual
        )
        summaries[[model_label]] <- data.frame(
          AnalysisSet = analysis_set, Model = model_label,
          Category = levels(droplevels(subset_data$Category)),
          N = as.integer(counts), Estimate = boot$estimate,
          CI95Low = boot$low, CI95High = boot$high,
          BootstrapReplicates = n_boot, Seed = seed + 100L,
          IdentityReferencePercent = design$centers[1],
          FlyCoverageReferencePercent = design$centers[2],
          HumanCoverageReferencePercent = design$centers[3]
        )
      }
    }
  }
  tests <- bind_rows(tests)
  tests$PAdjustedBH <- p.adjust(tests$PValue, method = "BH")
  tests$AdjustmentFamily <- paste(
    analysis_set, ": four tests (global and disease-only vs other, unadjusted and adjusted)"
  )
  list(
    tests = tests, summaries = bind_rows(summaries),
    coefficients = bind_rows(coefficients), diagnostics = bind_rows(diagnostics),
    gene_level_values = bind_rows(fitted_values)
  )
}

all_results <- run_set(values, "All measured orthologues", 9000L)
for (component in c("tests", "summaries", "coefficients", "diagnostics", "gene_level_values")) {
  table <- all_results[[component]]
  write.csv(table, summary_path(paste0("FigS8_", component, ".csv")), row.names = FALSE)
}
sources <- c(
  conservation_file, input_path("Assoc_Data.txt"), output_path("dmel_gene_reference.csv"),
  output_path("sexspecific_disease_ids.csv"), output_path("all_disease_ortholog_ids.csv"),
  output_path("all_ortholog_ids.csv"),
  file.path(project_root, "Sup_plotting_and_analyses", "figureS8_conservation_adjusted_sa.R"),
  file.path(project_root, "Sup_plotting_and_analyses", "conservation_sa_regression.R")
)
manifest <- data.frame(
  Input = substring(sources, nchar(project_root) + 2L),
  Bytes = file.info(sources)$size, MD5 = unname(tools::md5sum(sources))
)
write.csv(manifest, summary_path("FigS8_input_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), summary_path("FigS8_R_session.txt"))

plot_data <- all_results$summaries |>
  mutate(
    Category = factor(Category, levels = category_order),
    Model = factor(Model, levels = c("Unadjusted (same genes)", "Conservation adjusted"))
  )
offset <- position_dodge(width = 0.55)
figure <- ggplot(plot_data, aes(Model, Estimate, colour = Category, group = Category)) +
  geom_line(position = offset, linewidth = 0.45, alpha = 0.65) +
  geom_errorbar(aes(ymin = CI95Low, ymax = CI95High),
                position = offset, width = 0.12, linewidth = 0.65) +
  geom_point(position = offset, size = 2.8) +
  geom_text(aes(y = CI95High, label = paste0("n=", N)), position = offset,
            colour = "black", size = 3, vjust = -0.8) +
  scale_colour_manual(values = palette, breaks = category_order, drop = FALSE,
    labels = c("SSDOs", "Disease orthologues", "Orthologues")) +
  scale_x_discrete(labels = c(
    "Unadjusted (same genes)" = "Unadjusted\n(same genes)",
    "Conservation adjusted" = "Conservation\nadjusted"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(x = NULL, y = "SA effect", colour = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(colour = "black"),
    legend.position = "bottom", legend.text = element_text(size = 10),
    strip.background = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    axis.line = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1))
ggsave(file.path(project_root, "Sup_figures", "FigS8_conservation_adjusted_sa.png"),
       figure, width = 8.2, height = 5.7, dpi = 300, bg = "white")
message("Wrote FigS8, matched-universe estimates, clustered regression tests and input audits.")
