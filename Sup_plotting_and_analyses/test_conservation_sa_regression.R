#!/usr/bin/env Rscript
source("Sup_plotting_and_analyses/conservation_sa_regression.R")

set.seed(12)
g <- rep(seq_len(60), each = 2)
z <- rnorm(length(g))
category <- factor(rep(c("A", "B", "C"), each = 40))
y <- 0.7 * z + c(0, 0.4, 1)[category] +
  rep(rnorm(60, sd = 0.4), each = 2) + rnorm(length(g), sd = 0.2)
x <- model.matrix(~ z + category)
x0 <- model.matrix(~ z)
h <- matrix(0, 2, ncol(x))
h[1, 3] <- 1
h[2, 4] <- 1
fit <- fit_conservation_model(x, y, g)
expected <- lm(y ~ z + category)
stopifnot(isTRUE(all.equal(unname(coef(expected)), unname(fit$beta), tolerance = 1e-10)))

scores <- rowsum(x * residuals(expected), g)
inverse <- solve(crossprod(x))
v <- 60 / 59 * 119 / 116 * inverse %*% crossprod(scores) %*% inverse
stopifnot(isTRUE(all.equal(unname(v), unname(fit$covariance), tolerance = 1e-10)))
test <- wild_cluster_test(x, x0, y, g, h, n_resamples = 199L, seed = 7L)
repeat_test <- wild_cluster_test(x, x0, y, g, h, n_resamples = 199L, seed = 7L)
stopifnot(
  identical(test$null_statistics, repeat_test$null_statistics),
  test$p_value == (sum(test$null_statistics >= test$statistic - 1e-12) + 1) / 200,
  test$p_value <= 0.05
)
set.seed(7L)
signs <- sample(c(-1, 1), 60, replace = TRUE)
null_fit <- lm(y ~ z)
y_star <- fitted(null_fit) + residuals(null_fit) * signs[g]
explicit_star <- fit_conservation_model(x, y_star, g)
stopifnot(isTRUE(all.equal(
  test$null_statistics[1],
  wald_statistic(explicit_star$beta, explicit_star$covariance, h),
  tolerance = 1e-10
)))

# Covariate adjustment is invariant to units and affine shifts of a covariate.
x_scaled <- model.matrix(~ I(10 + 3 * z) + category)
scaled <- wild_cluster_test(x_scaled, model.matrix(~ I(10 + 3 * z)), y, g, h,
                           n_resamples = 199L, seed = 7L)
stopifnot(isTRUE(all.equal(test$statistic, scaled$statistic, tolerance = 1e-9)),
          test$p_value == scaled$p_value)

prediction <- model.matrix(~ z + category,
  data.frame(z = 0, category = factor(c("A", "B", "C"), levels = levels(category))))
boot <- cluster_bootstrap_predictions(x, y, g, prediction, n_boot = 199L, seed = 4L)
stopifnot(
  isTRUE(all.equal(boot$estimate, drop(prediction %*% coef(expected)))),
  all(boot$low < boot$high), all(is.finite(boot$replicates))
)
# Check the optimized cluster bootstrap against an explicit replicated-row fit.
set.seed(4L)
sampled <- sample.int(60, 60, replace = TRUE)
indices <- unlist(lapply(sampled, function(i) which(g == i)))
explicit <- lm.fit(x[indices, , drop = FALSE], y[indices])
stopifnot(isTRUE(all.equal(
  boot$replicates[1, ], drop(prediction %*% explicit$coefficients),
  check.attributes = FALSE, tolerance = 1e-9
)))
bad <- try(fit_conservation_model(cbind(x, x[, 2]), y, g), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
cat("Conservation-adjusted regression tests passed.\n")
