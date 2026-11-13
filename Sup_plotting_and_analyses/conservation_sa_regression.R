# Small fixed-design regressions with human-partner-clustered inference.
fit_conservation_model <- function(x, y, cluster) {
  if (!is.matrix(x) || !is.numeric(x) || length(y) != nrow(x) ||
      length(cluster) != length(y) || anyNA(cluster) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Expected a finite design/response and one human-partner cluster per row.")
  }
  n <- nrow(x)
  p <- ncol(x)
  cluster <- match(cluster, unique(cluster))
  g <- max(cluster)
  if (qr(x)$rank != p || n <= p || g <= p) {
    stop("Insufficient independent clusters or rank-deficient regression design.")
  }
  inverse <- solve(crossprod(x))
  beta <- drop(inverse %*% crossprod(x, y))
  fitted <- drop(x %*% beta)
  residual <- y - fitted
  score <- rowsum(x * residual, cluster, reorder = FALSE)
  correction <- g / (g - 1) * (n - 1) / (n - p)
  covariance <- correction * inverse %*% crossprod(score) %*% inverse
  list(
    beta = beta, fitted = fitted, residual = residual, covariance = covariance,
    inverse = inverse, cluster = cluster, n = n, p = p, g = g,
    correction = correction
  )
}

wald_statistic <- function(beta, covariance, hypothesis) {
  effect <- drop(hypothesis %*% beta)
  variance <- hypothesis %*% covariance %*% t(hypothesis)
  if (any(!is.finite(variance)) || qr(variance)$rank != nrow(hypothesis)) {
    stop("Singular hypothesis covariance.")
  }
  drop(crossprod(effect, solve(variance, effect))) / nrow(hypothesis)
}

wild_cluster_test <- function(x, reduced, y, cluster, hypothesis,
                              n_resamples = 10000L, seed = 1L) {
  full_fit <- fit_conservation_model(x, y, cluster)
  reduced_fit <- fit_conservation_model(reduced, y, cluster)
  if (qr(cbind(x, reduced))$rank != ncol(x)) {
    stop("The null design must be nested in the full design.")
  }
  if (max(abs(hypothesis %*% full_fit$inverse %*%
      crossprod(x, reduced_fit$fitted))) > 1e-8) {
    stop("The null design does not impose the specified hypothesis.")
  }
  observed <- wald_statistic(full_fit$beta, full_fit$covariance, hypothesis)
  set.seed(seed)
  simulated <- numeric(n_resamples)
  for (b in seq_len(n_resamples)) {
    signs <- sample(c(-1, 1), full_fit$g, replace = TRUE)
    y_star <- reduced_fit$fitted +
      reduced_fit$residual * signs[full_fit$cluster]
    beta_star <- drop(full_fit$inverse %*% crossprod(x, y_star))
    residual_star <- y_star - drop(x %*% beta_star)
    scores <- rowsum(x * residual_star, full_fit$cluster, reorder = FALSE)
    covariance_star <- full_fit$correction * full_fit$inverse %*%
      crossprod(scores) %*% full_fit$inverse
    simulated[b] <- wald_statistic(beta_star, covariance_star, hypothesis)
  }
  list(
    fit = full_fit, statistic = observed, null_statistics = simulated,
    p_value = (sum(simulated >= observed - 1e-12) + 1) / (n_resamples + 1),
    seed = seed, resamples = n_resamples
  )
}

cluster_bootstrap_predictions <- function(x, y, cluster, prediction,
                                          n_boot = 5000L, seed = 1L) {
  fit <- fit_conservation_model(x, y, cluster)
  p <- fit$p
  group <- fit$cluster
  # Cluster sufficient statistics avoid refitting thousands of replicated gene tables.
  xy <- rowsum(x * y, group, reorder = FALSE)
  xx <- matrix(0, fit$g, p * p)
  for (j in seq_len(p)) {
    for (k in seq_len(p)) {
      xx[, j + (k - 1L) * p] <- drop(rowsum(x[, j] * x[, k], group, reorder = FALSE))
    }
  }
  set.seed(seed)
  draws <- matrix(NA_real_, n_boot, nrow(prediction))
  for (b in seq_len(n_boot)) {
    weights <- tabulate(sample.int(fit$g, fit$g, replace = TRUE), nbins = fit$g)
    xx_star <- matrix(drop(crossprod(weights, xx)), p, p)
    xy_star <- drop(crossprod(weights, xy))
    if (qr(xx_star)$rank != p) {
      stop("Cluster bootstrap lost a category or produced a singular design.")
    }
    draws[b, ] <- drop(prediction %*% solve(xx_star, xy_star))
  }
  list(
    estimate = drop(prediction %*% fit$beta),
    low = apply(draws, 2, quantile, probs = 0.025, names = FALSE),
    high = apply(draws, 2, quantile, probs = 0.975, names = FALSE),
    replicates = draws, seed = seed, n_boot = n_boot
  )
}
