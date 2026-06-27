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