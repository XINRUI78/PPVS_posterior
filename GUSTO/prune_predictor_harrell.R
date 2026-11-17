prune_predictor_harrell <- function(x, eta_full, r2_threshold = 0.99) {
  
  x <- as.data.frame(x)
  eta_full <- as.numeric(eta_full)
  predictors <- names(x)
  p <- ncol(x)
  
  empty_b <- c("Intercept" = mean(eta_full), setNames(rep(0, p), predictors))
  empty_kept <- setNames(integer(p), predictors)
  
  # Fit auxiliary OLS using estimable predictors
  dat <- data.frame(eta_full = eta_full, x[, predictors, drop = FALSE])
  fm <- reformulate(predictors, response = "eta_full")
  
  a <- rms::ols(fm, data = dat, sigma = 1, x = TRUE, y = TRUE)
  s <- rms::fastbw(a, aics = 10000)
  
  # Extract auxiliary OLS coefficients
  X <- model.matrix(fm, data = dat)
  cn <- colnames(X)
  cn[1] <- "Intercept"
  full_betas <- coef(a)[cn]
  
  # Extract fastbw coefficient history
  B <- as.matrix(s$Coefficients)
  
  # Include full auxiliary model
  B <- rbind(full_betas, B)
  
  # Approximate linear predictors
  eta_path <- X %*% t(B)
  
  # R-squared for each step
  r2 <- vapply(seq_len(ncol(eta_path)), function(j) {
    z <- eta_path[, j]
    if (sd(z) < 1e-12) return(0)
    cor(z, eta_full)^2
  }, numeric(1))
  
  path <- data.frame(step = seq_along(r2) - 1L, r2 = r2)
  
  # Select the smallest model with R-squared >= threshold
  best <- max(which(r2 >= r2_threshold))
  
  b <- B[best, ]
  
  # Restore coefficients for all original predictors
  b_complete <- empty_b
  b_complete[names(b)] <- b
  
  selected <- names(which(abs(b_complete[predictors]) > 1e-12))
  kept <- as.integer(paste0("X", 1:30) %in% selected)
  
  list(coefficients = b_complete, kept = kept,
       r2 = unname(r2[best]), path = path)
}

# ============================================================
# Validation and fixed-patient predictions
# ============================================================

predict_harrell_stepdown <- function(result, newdata) {
  b <- result$coefficients
  X <- cbind("Intercept" = 1, as.matrix(newdata[,  names(b)[-1], drop = FALSE]))
  plogis(as.vector(X %*% b))
}
