library(tidyverse)
library(rstanarm)
library(projpred)
library(loo)
library(glmnet)
library(foreach)
library(doParallel)
library(MASS)

run_hs_p <- function(
    i,
    ndev,
    n.para,
    beta0,
    beta,
    fixed_patients,
    n.true = NULL,
    nterms_max = 30,
    prediction_draws = 400,
    projection_draws = 400
) {
  
  set.seed(i)
  SEED <- as.integer(i)
  
  # Store full draws only for simulations 1 to 10
  save_full_draws <- i <= 10
  
  # ==========================================================
  # Extract saved patient information
  # ==========================================================
  
  patient_y <- fixed_patients$y
  p_fixed_true <- fixed_patients$p_true
  
  predictor_names <- paste0(
    "X",
    seq_len(n.para)
  )
  
  fixed_patients_x <- fixed_patients[,predictor_names,drop = FALSE]
  
  # ==========================================================
  # Patient names
  # ==========================================================
    patient_names <- c(
      "low1",
      "low2",
      "medium1",
      "medium2",
      "high1",
      "high2"
    )
 
    rownames(fixed_patients) <- patient_names 
  
  # ==========================================================
  # Generate development dataset
  # ==========================================================
  
  data.dev <- generate_ss(
    n = ndev,
    n.para = n.para,
    beta0 = beta0,
    beta = beta,
    n.true = n.true
  )
  
  all_vars <- colnames(data.dev)[-1]
  
  # Put fixed-patient predictors in the same order as data.dev
  fixed_patients_x <- fixed_patients_x[
    ,
    all_vars,
    drop = FALSE
  ]
  
  # ==========================================================
  # Fit horseshoe reference model
  # ==========================================================
  
  ref_formula <- as.formula(
    paste("y", paste(all_vars, collapse = " + "), sep = " ~ ")
  )
  
  ref_fit <- stan_glm(
    formula = ref_formula,
    data = data.dev,
    family = binomial(link = "logit"),
    prior = hs(),
    prior_intercept = normal(0, 1),
    QR = TRUE,
    seed = SEED + 2000,
    adapt_delta = 0.99,
    iter = 4000,
    chains = 4,
    cores = 4,
    refresh = 0
  )
  
  # Keep information required by projpred
  ref_fit$call$formula <- ref_formula
  ref_fit$call$data <- data.dev
  ref_fit$call$seed <- SEED + 2000
  
  # ==========================================================
  # Projection-predictive variable selection
  # ==========================================================
  
  vs <- cv_varsel(
    object = ref_fit,
    method = "forward",
    cv_method = "kfold",
    K = 5,
    validate_search = TRUE,
    nterms_max = min(
      nterms_max,
      length(all_vars)
    ),
    seed = SEED + 3000
  )
  
  # Suggested model size
  size_suggested <- suggest_size(
    vs
  )
  
  # Best ELPD model size
  get_best_size <- function(vs) {
    perf <- performances(vs)$submodels
    perf$size[which.max(replace(perf$elpd, is.na(perf$elpd), -Inf))]
  }
  
  size_best <- get_best_size(vs)
  
  # ==========================================================
  # Project onto suggested and best submodels
  # ==========================================================
  
  proj_suggest <- project(
    object = vs,
    nterms = size_suggested,
    ns = projection_draws,
    seed = SEED + 40000
  )
  
  proj_best <- project(
    object = vs,
    nterms = size_best,
    ns = projection_draws,
    seed = SEED + 40001
  )
  
  # ==========================================================
  # Posterior probability draws for the six fixed patients
  # ==========================================================
  
  # ----------------------------------------------------------
  # Horseshoe reference model
  # ----------------------------------------------------------
  
  eta_ref_patients <- posterior_linpred(
    object = ref_fit,
    newdata = fixed_patients_x,
    transform = FALSE,
    draws = prediction_draws
  )
  
  eta_ref_patients <- orient_draw_matrix(
    draw_matrix = eta_ref_patients,
    n_observations = nrow(fixed_patients_x)
  )
  
  p_ref_patients <- plogis(
    eta_ref_patients
  )
  
  colnames(p_ref_patients) <- patient_names
  
  # ----------------------------------------------------------
  # Suggested projected model
  # ----------------------------------------------------------
  
  eta_suggest_patients <- proj_linpred(
    object = proj_suggest,
    newdata = fixed_patients_x,
    integrated = FALSE,
    transform = FALSE
  )$pred
  
  eta_suggest_patients <- orient_draw_matrix(
    draw_matrix = eta_suggest_patients,
    n_observations = nrow(fixed_patients_x)
  )
  
  p_suggest_patients <- plogis(
    eta_suggest_patients
  )
  
  colnames(p_suggest_patients) <- patient_names
  
  # ----------------------------------------------------------
  # Best projected model
  # ----------------------------------------------------------
  
  eta_best_patients <- proj_linpred(
    object = proj_best,
    newdata = fixed_patients_x,
    integrated = FALSE,
    transform = FALSE
  )$pred
  
  eta_best_patients <- orient_draw_matrix(
    draw_matrix = eta_best_patients,
    n_observations = nrow(fixed_patients_x)
  )
  
  p_best_patients <- plogis(
    eta_best_patients
  )
  
  colnames(p_best_patients) <- patient_names
  
  # ==========================================================
  # Task 1
  #
  # Store posterior mean probability for each patient in every
  # simulation.
  # ==========================================================
  
  patient_mean_predictions <- bind_rows(
    
    tibble(
      simulation = i,
      ndev = ndev,
      method = "hs-ref",
      patient = patient_names,
      observed_y = as.integer(patient_y),
      true_probability = as.numeric(p_fixed_true),
      predicted_probability = colMeans(
        p_ref_patients
      )
    ),
    
    tibble(
      simulation = i,
      ndev = ndev,
      method = "hs-1se",
      patient = patient_names,
      observed_y = as.integer(patient_y),
      true_probability = as.numeric(p_fixed_true),
      predicted_probability = colMeans(
        p_suggest_patients
      )
    ),
    
    tibble(
      simulation = i,
      ndev = ndev,
      method = "hs-best",
      patient = patient_names,
      observed_y = as.integer(patient_y),
      true_probability = as.numeric(p_fixed_true),
      predicted_probability = colMeans(
        p_best_patients
      )
    )
  )
  
  # ==========================================================
  # Task 2
  #
  # Store full posterior probability and predictive draws only
  # for simulations 1 to 10.
  # ==========================================================
  
  patient_draws <- NULL
  
  if (save_full_draws) {
    
    # --------------------------------------------------------
    # Reference-model posterior predictive draws
    # --------------------------------------------------------
    
    yrep_ref_patients <- posterior_predict(
      object = ref_fit,
      newdata = fixed_patients_x,
      draws = prediction_draws
    )
    
    yrep_ref_patients <- orient_draw_matrix(
      draw_matrix = yrep_ref_patients,
      n_observations = nrow(fixed_patients_x)
    )
    
    colnames(yrep_ref_patients) <- patient_names
    
    # --------------------------------------------------------
    # Suggested-model posterior predictive draws
    # --------------------------------------------------------
    
    set.seed(SEED + 50000)
    
    yrep_suggest_patients <- matrix(
      rbinom(
        n = length(p_suggest_patients),
        size = 1,
        prob = as.vector(p_suggest_patients)
      ),
      nrow = nrow(p_suggest_patients),
      ncol = ncol(p_suggest_patients),
      dimnames = list(
        NULL,
        patient_names
      )
    )
    
    # --------------------------------------------------------
    # Best-model posterior predictive draws
    # --------------------------------------------------------
    
    set.seed(SEED + 50001)
    
    yrep_best_patients <- matrix(
      rbinom(
        n = length(p_best_patients),
        size = 1,
        prob = as.vector(p_best_patients)
      ),
      nrow = nrow(p_best_patients),
      ncol = ncol(p_best_patients),
      dimnames = list(
        NULL,
        patient_names
      )
    )
    
    patient_draws <- list(
      simulation = i,
      ndev = ndev,
      patient_names = patient_names,
      observed_y = patient_y,
      true_probability = p_fixed_true,
      
      probability_draws = list(
        hs_ref = p_ref_patients,
        hs_1se = p_suggest_patients,
        hs_best = p_best_patients
      ),
      
      predictive_draws = list(
        hs_ref = yrep_ref_patients,
        hs_1se = yrep_suggest_patients,
        hs_best = yrep_best_patients
      )
    )
  }
  
  # ==========================================================
  # Return results
  # ==========================================================
  
  list(
    patient_mean_predictions = patient_mean_predictions,
    patient_draws = patient_draws
  )
}
