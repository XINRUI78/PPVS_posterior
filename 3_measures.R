library(speedglm)
measures <- function(yval, p_val) {
  eta_val <- log(p_val/(1 - p_val))
  # Calibration slope
  fitcal <- speedglm(yval ~ eta_val, family = binomial())
  cal_slope <- as.vector(coef(fitcal)[2])
  # Calibration in the large
  off <- speedglm(yval ~ 1, offset = eta_val, family = binomial())
  cal_large <- as.vector(coef(off))
  # AUC
  cstat <- roc(response = yval, predictor = as.vector(p_val), levels = c(0, 1), direction = "<")
  auc <- as.vector(cstat$auc)
  # Root mean square prediction error (RMSPE)
  rmspe <- sqrt(mean((p_val - yval)^2))
  return(c(cal_slope, cal_large, auc, rmspe))
}

# Performance measures over posterior draws
# Default: return only mean over draws
# If return_all = TRUE: return both mean and all draw-level measures
measures_by_draw <- function(yval,
                             p_draws,
                             return_all = FALSE) {
  
  # p_draws: draws × observations
  
  out <- t(apply(p_draws, 1, function(p) {
    measures(yval, p)
  }))
  
  colnames(out) <- c(
    "calibration_slope",
    "calibration_in_the_large",
    "auc",
    "rmspe"
  )
  
  mean_measures <- colMeans(out, na.rm = TRUE)
  
  if (!return_all) {
    return(mean_measures)
  }
  
  return(list(
    mean = mean_measures,
    draws = as.data.frame(out)
  ))
}