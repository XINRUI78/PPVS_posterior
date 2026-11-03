prune_lasso_predictor <- function(x, eta_lasso, r2_threshold = 0.99) {
  x <- as.data.frame(x)
  eta_lasso <- as.numeric(eta_lasso)
  predictors <- names(x)
  dat <- data.frame(eta_lasso = eta_lasso, x, check.names = FALSE)
  fit <- function(terms) {
    fm <- if (length(terms)) reformulate(terms, response = "eta_lasso")
          else eta_lasso ~ 1
    lm(fm, data = dat)
  }
  # For a constant target, the intercept-only model reproduces it exactly.
  if (stats::var(eta_lasso) < 1e-20) {
    result <- fit(character(0))
    return(list(coefficients = coef(result),
                kept = setNames(integer(length(predictors)), predictors),
                r2 = 1))
  }
  remaining <- predictors
  while (length(remaining) > 0L) {
    candidates <- lapply(remaining, function(drop) {
      subset <- setdiff(remaining, drop)
      mod <- fit(subset)
      r2 <- summary(mod)$r.squared
      if (length(subset) == 0L) r2 <- 0
      if (!is.finite(r2)) r2 <- -Inf
      list(drop = drop, r2 = r2)
    })
    scores <- vapply(candidates, `[[`, numeric(1), "r2")
    j <- which.max(scores)
    if (!length(j) || scores[j] < r2_threshold) break
    remaining <- setdiff(remaining, candidates[[j]]$drop)
  }
  final <- fit(remaining)
  b <- coef(final)
  if (anyNA(b)) stop("Aliased coefficients in pruned model; check collinearity")
  kept <- setNames(as.integer(predictors %in% remaining), predictors)
  list(coefficients = b, kept = kept,
       r2 = if (length(remaining)) summary(final)$r.squared else 0)
}
