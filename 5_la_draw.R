library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

run_la <- function(
    i,
    ndev,
    nval,
    n.para,
    beta0,
    beta,
    n.true = NULL,
    nterms_max = 30,
    prediction_draws = 400,
    projection_draws = 400
) {
  
  set.seed(i)
  SEED <- as.integer(i)
  
  # Detailed calculations are performed only for:
  # 1. the first 10 simulations; and
  # 2. ndev >= 1500
  run_draw_details <- i <= 10 && ndev >= 1500
  
  # ----------------------------------------------------------
  # Generate data
  # ----------------------------------------------------------
  
  data.dev <- generate_ss(
    n = ndev,
    n.para = n.para,
    beta0 = beta0,
    beta = beta,
    n.true = n.true
  )
  
  data.val <- generate_ss(
    n = nval,
    n.para = n.para,
    beta0 = beta0,
    beta = beta,
    n.true = n.true
  )
  
  xval <- data.val[, -1]
  yval <- as.numeric(data.val[, 1])
  
  # True validation probabilities

  eta_true <- as.numeric(
    beta0 +
      as.matrix(xval) %*% beta
  )
  
  p_true <- plogis(eta_true)
  
  # Fit Laplace-prior reference model

  all_vars <- colnames(data.dev)[-1]
  
  ref_formula <- as.formula(
    paste("y", paste(all_vars, collapse = " + "), sep = " ~ ")
  )
  
  ref_fit <- stan_glm(
    formula = ref_formula,
    data = data.dev,
    family = binomial(link = "logit"),
    prior = laplace(location = 0, scale = 1),
    prior_intercept = normal(0, 1),
    QR = TRUE,
    seed = SEED + 2000,
    adapt_delta = 0.99,
    iter = 4000,
    cores = 4,
    chains = 4
  )
  
  ref_fit$call$formula <- ref_formula
  ref_fit$call$data <- data.dev
  ref_fit$call$seed <- SEED + 2000
  
  # ----------------------------------------------------------
  # Reference-model linear-predictor draws
  # ----------------------------------------------------------
  
  eta_ref_draws <- posterior_linpred(
    object = ref_fit,
    newdata = xval,
    transform = FALSE,
    draws = prediction_draws
  )
  
  eta_ref_draws <- orient_draw_matrix(
    draw_matrix = eta_ref_draws,
    n_observations = length(yval)
  )
  
  ref_predictions <- summarise_prediction_draws(
    eta_draws = eta_ref_draws
  )
  
  # Projection predictive selection

  vs <- cv_varsel(
    object = ref_fit,
    method = "forward",
    cv_method = "kfold",
    K = 5,
    validate_search = TRUE,
    nterms_max = min(
      nterms_max,
      n.para
    ),
    seed = SEED + 3000
  )
  
  # Selected-variable indicators
  
  get_selected <- function(vs, k) {
    sel_vars <- ranking(vs)$fulldata[1:k]
    as.integer(all_vars %in% sel_vars)
  }
 
  # suggest model size
  size_suggested <- suggest_size(vs)
  
  # Best ELPD model size
  get_best_size <- function(vs) {
    perf <- performances(vs)$submodels
    perf$size[which.max(replace(perf$elpd, is.na(perf$elpd), -Inf))]
  }
  
  size_best <- get_best_size(vs)
  
  # Project onto suggested and best models
  proj_suggest <- project(
    vs,
    nterms = size_suggested,
    ns = projection_draws,
    seed = SEED + 40000
  )
  
  proj_best <- project(
    vs,
    nterms = size_best,
    ns = projection_draws,
    seed = SEED + 40001
  )
  
  # Suggested-model linear-predictor draws
  
  eta_suggest_draws <- proj_linpred(
    proj_suggest,
    newdata = xval,
    integrated = FALSE,
    transform = FALSE
  )$pred
  
  eta_suggest_draws <- orient_draw_matrix(
    draw_matrix = eta_suggest_draws,
    n_observations = length(yval)
  )
  
  suggest_predictions <- summarise_prediction_draws(
    eta_draws = eta_suggest_draws
  )
  
  # ----------------------------------------------------------
  # Best-model linear-predictor draws
  # ----------------------------------------------------------
  
  eta_best_draws <- proj_linpred(
    proj_best,
    newdata = xval,
    integrated = FALSE,
    transform = FALSE
  )$pred
  
  eta_best_draws <- orient_draw_matrix(
    draw_matrix = eta_best_draws,
    n_observations = length(yval)
  )
  
  best_predictions <- summarise_prediction_draws(
    eta_draws = eta_best_draws
  )
  
  # ----------------------------------------------------------
  # Main summary
  #
  # 1 ndev column
  # 16 performance columns
  # n.para variable-selection columns
  # ----------------------------------------------------------
  
  ref_row <- c(
    ndev, "la-ref",
    measures(yval, ref_predictions$mean_beta),
    measures(yval, ref_predictions$mode_beta),
    measures(yval, ref_predictions$mean_probability),
    measures(yval, ref_predictions$mode_probability),
    rep(1, n.para)
  )
  
  suggest_row <- c(
    ndev, "la-1se",
    measures(yval, suggest_predictions$mean_beta),
    measures(yval, suggest_predictions$mode_beta),
    measures(yval, suggest_predictions$mean_probability),
    measures(yval, suggest_predictions$mode_probability),
    get_selected(vs, size_suggested)
  )
  
  best_row <- c(
    ndev, "la-best",
    measures(yval, best_predictions$mean_beta),
    measures(yval, best_predictions$mode_beta),
    measures(yval, best_predictions$mean_probability),
    measures(yval, best_predictions$mode_probability),
    get_selected(vs, size_best)
  )
  
  summary_result <- rbind(
    ref_row,
    suggest_row,
    best_row
  )
  
  measure_names <- c(
    "calibration_slope",
    "calibration_in_the_large",
    "auc",
    "rmspe"
  )
  
  colnames(summary_result) <- c(
    "ndev", "method",
    paste0(measure_names, "_mean_beta"),
    paste0(measure_names, "_mode_beta"),
    paste0(measure_names, "_mean_probability"),
    paste0(measure_names, "_mode_probability"),
    paste0("varsel", seq_len(n.para))
  )
  
  summary_result <- as.data.frame(
    summary_result,
    check.names = FALSE
  )
  
  # Stop here when detailed calculations are not required
  
  if (!run_draw_details) {
    return(list(summary = summary_result))
  }
  
  # Everything below is completely skipped when ndev < 1500
  # or when i > 10.
  
  # Draw-specific performance measures
  
  ref_draw_measures <- measures_by_draw(
    yval = yval,
    p_draws = ref_predictions$p_draws
  ) %>%
    mutate(
      method = "la-ref",
      simulation = i,
      ndev = ndev
    )
  
  suggest_draw_measures <- measures_by_draw(
    yval = yval,
    p_draws = suggest_predictions$p_draws
  ) %>%
    mutate(
      method = "la-1se",
      simulation = i,
      ndev = ndev
    )
  
  best_draw_measures <- measures_by_draw(
    yval = yval,
    p_draws = best_predictions$p_draws
  ) %>%
    mutate(
      method = "la-best",
      simulation = i,
      ndev = ndev
    )
  
  draw_measures_all <- bind_rows(
    ref_draw_measures,
    suggest_draw_measures,
    best_draw_measures
  ) 
  
  # ----------------------------------------------------------
  # Six patients selected using true risk
  # ----------------------------------------------------------
 
  selected_ids <- select_patients(
    p_true = p_true
  )
  
  patient_names <- c(
    "low1",
    "low2",
    "medium1",
    "medium2",
    "high1",
    "high2"
  )
  
  names(selected_ids) <- patient_names
  
  patient_prob_draws <- list(
    patient_ids = selected_ids,
    yval = yval[selected_ids],
    ptrue = p_true[selected_ids],
    ref_la = make_patient_draws(ref_predictions$p_draws, selected_ids, length(yval)),
    la_1se = make_patient_draws(suggest_predictions$p_draws, selected_ids, length(yval)),
    la_best = make_patient_draws(best_predictions$p_draws, selected_ids, length(yval))
  )
  
  list(
    summary = summary_result,
    draw_measures_all = draw_measures_all,
    patient_prob_draws = patient_prob_draws
  )
}
  
 
