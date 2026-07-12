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

# ============================================================

orient_draw_matrix <- function(draw_matrix, n_observations) {
  
  draw_matrix <- as.matrix(draw_matrix)
  
  if (ncol(draw_matrix) == n_observations) {
    return(draw_matrix)
  }
  
  if (nrow(draw_matrix) == n_observations) {
    return(t(draw_matrix))
  }
  
  stop(
    "Neither dimension of the prediction matrix matches ",
    "the number of validation observations."
  )
}

# ============================================================
# Four summaries of posterior predictions
# mean_beta:
#   inverse-logit of mean linear predictor
# mode_beta:
#   inverse-logit of modal linear predictor
# mean_probability:
#   mean probability across posterior draws
# mode_probability:
#   modal probability across posterior draws

summarise_prediction_draws <- function(eta_draws) {
  
  eta_draws <- as.matrix(eta_draws)
  
  p_draws <- plogis(eta_draws)
  
  eta_mean <- colMeans(
    eta_draws,
    na.rm = TRUE
  )
  
  eta_mode <- apply(
    eta_draws,
    2,
    estimate_mode
  )
  
  p_mean <- colMeans(
    p_draws,
    na.rm = TRUE
  )
  
  p_mode <- apply(
    p_draws,
    2,
    estimate_mode
  )
  
  list(
    p_draws = p_draws,
    mean_beta = plogis(eta_mean),
    mode_beta = plogis(eta_mode),
    mean_probability = p_mean,
    mode_probability = p_mode
  )
}


# ============================================================
# 6. Measures for every posterior draw
# ============================================================

measures_by_draw <- function(yval, p_draws) {
  
  p_draws <- orient_draw_matrix(
    draw_matrix = p_draws,
    n_observations = length(yval)
  )
  
  out <- t(
    apply(
      p_draws,
      1,
      function(p) {
        measures(yval, p)
      }
    )
  )
  
  out <- as.data.frame(out)
  
  colnames(out) <- c(
    "calibration_slope",
    "calibration_in_the_large",
    "auc",
    "rmspe"
  )
  
  out$draw <- seq_len(nrow(out))
  
  out
}

