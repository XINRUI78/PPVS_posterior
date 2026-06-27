library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

run_la <- function(i, ndev, nval, n.para, beta0, beta,
                   nterms_max = 30,
                   ns = 2000,
                   save_draw_details = i <= 10) {
  
  set.seed(i)
  SEED <- as.integer(i)
  
  data.dev <- generate_ss(ndev, n.para, beta0, beta)
  data.val <- generate_ss(nval, n.para, beta0, beta)
  
  xval <- data.val[, -1]
  yval <- data.val[, 1]
  
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
  
  # Reference model posterior probability draws
  p_ref_draws <- posterior_epred(ref_fit, newdata = xval)
  
  # Mean probability for original summary performance
  p_ref <- colMeans(p_ref_draws)
  
  ref_draw_results <- measures_by_draw(
    yval = yval,
    p_draws = p_ref_draws,
    return_all = save_draw_details
  )
  
  ref_draw_mean <- if (save_draw_details) {
    ref_draw_results$mean
  } else {
    ref_draw_results
  }
  
  selected_patients <- select_patients(p_ref)
  
  ref_row <- c(
    0.3, 0.8, ndev, "ref-la",
    measures(yval, p_ref), ref_draw_mean,
    rep(1, n.para),
    NA
  )
  
  vs <- cv_varsel(
    ref_fit,
    method = "forward",
    cv_method = "kfold",
    K = 5,
    validate_search = TRUE,
    nterms_max = min(nterms_max, n.para),
    seed = SEED + 3000
  )
  
  get_selected <- function(vs, k) {
    sel_vars <- ranking(vs)$fulldata[1:k]
    as.integer(all_vars %in% sel_vars)
  }
  
  size_suggested <- suggest_size(vs)
  
  get_best_size <- function(vs) {
    perf <- performances(vs)$submodels
    perf$size[which.max(replace(perf$elpd, is.na(perf$elpd), -Inf))]
  }
  
  size_best <- get_best_size(vs)
  
  proj_suggest <- project(
    vs,
    nterms = size_suggested,
    seed = SEED + 40000,
    ns = ns
  )
  
  proj_best <- project(
    vs,
    nterms = size_best,
    seed = SEED + 40001,
    ns = ns
  )
  
  p_suggest_draws <- proj_linpred(
    proj_suggest,
    newdata = xval,
    integrated = FALSE,
    transform = TRUE
  )$pred
  
  p_best_draws <- proj_linpred(
    proj_best,
    newdata = xval,
    integrated = FALSE,
    transform = TRUE
  )$pred
  
  if (ncol(p_suggest_draws) != length(yval)) {
    p_suggest_draws <- t(p_suggest_draws)
  }
  
  if (ncol(p_best_draws) != length(yval)) {
    p_best_draws <- t(p_best_draws)
  }
  
  p_suggest <- colMeans(p_suggest_draws)
  p_best <- colMeans(p_best_draws)
  
  suggest_draw_results <- measures_by_draw(
    yval = yval,
    p_draws = p_suggest_draws,
    return_all = save_draw_details
  )
  
  best_draw_results <- measures_by_draw(
    yval = yval,
    p_draws = p_best_draws,
    return_all = save_draw_details
  )
  
  suggest_draw_mean <- if (save_draw_details) {
    suggest_draw_results$mean
  } else {
    suggest_draw_results
  }
  
  best_draw_mean <- if (save_draw_details) {
    best_draw_results$mean
  } else {
    best_draw_results
  }
  
  suggest_row <- c(
    0.3, 0.8, ndev, "la-1se",
    measures(yval, p_suggest),
    suggest_draw_mean,
    get_selected(vs, size_suggested),
    NA
  )
  
  best_row <- c(
    0.3, 0.8, ndev, "la-best",
    measures(yval, p_best),
    best_draw_mean,
    get_selected(vs, size_best),
    NA
  )
  
  summary_result <- rbind(
    ref_row,
    suggest_row,
    best_row
  )
  
  colnames(summary_result) <- c(
    "prevalence",
    "anticipated c-stat",
    "ndev",
    "method",
    "calibration slope",
    "calibration in the large",
    "auc",
    "rmspe",
    "mean draw calibration slope",
    "mean draw calibration in the large",
    "mean draw auc",
    "mean draw rmspe",
    paste0("varsel", 1:n.para),
    "option"
  )
    
  if (save_draw_details) {
    
    selected_ids <- select_patients(p_ref)
    
    draw_measures_all <- bind_rows(
      ref_draw_results$draws %>%
        mutate(draw = row_number(), method = "ref-la", sim = i, ndev = ndev),
      
      suggest_draw_results$draws %>%
        mutate(draw = row_number(), method = "la-1se", sim = i, ndev = ndev),
      
      best_draw_results$draws %>%
        mutate(draw = row_number(), method = "la-best", sim = i, ndev = ndev)
    )
    
    patient_prob_draws <- list(
      patient_ids = selected_ids,
      yval = yval[selected_ids],
      ref_la = make_patient_draws(p_ref_draws, selected_ids, yval),
      la_1se = make_patient_draws(p_suggest_draws, selected_ids, yval),
      la_best = make_patient_draws(p_best_draws, selected_ids, yval)
    )
    
  } else {
    
    draw_measures_all <- NULL
    patient_prob_draws <- NULL
  }
  
return(list(
  summary = summary_result,
  draw_measures_all = draw_measures_all,
  patient_prob_draws = patient_prob_draws
))
}
  
 