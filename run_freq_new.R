library(glmnet)
library(rms)

# ============================================================
# Harrell step-down approximation
# ============================================================

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
# ============================================================
# Run one simulation replication
# ============================================================

run_freq_new <- function(i, scenario_id, ndev, n.para, beta0, beta, fixed_patients,
                         data_dir = "simulation_data") {
  
  # Read the existing simulation dataset
  data_file <- simulation_file(scenario_id, ndev, i, data_dir)
  sim_data <- readRDS(data_file)
  data.dev <- sim_data$data.dev
  data.val <- sim_data$data.val
  
  predictor_names <- paste0("X", seq_len(n.para))
  x <- data.dev[, predictor_names, drop = FALSE]
  y <- as.numeric(data.dev$y)
  xval <- data.val[, predictor_names, drop = FALSE]
  yval <- as.numeric(data.val$y)
  fixed_x <- fixed_patients[, predictor_names, drop = FALSE]
  
  x_mat <- as.matrix(x)
  xval_mat <- as.matrix(xval)
  fixed_mat <- as.matrix(fixed_x)
  p_true_val <- plogis(as.numeric(beta0 + xval_mat %*% as.numeric(beta)))
  p_true_pat <- as.numeric(fixed_patients$p_true)
  
  prev <- 0.3
  auc <- 0.8
  methods <- c("MLE", "UVS-5%", "UVS-15%", "BE-5%", "BE-15%", "LS-min",
               "LS-mod", "ALASSO", "SD-MLE", "SD-LS-min", "SD-LS-mod")
  method_result <- matrix(NA_character_, nrow = length(methods), ncol = 10 + n.para)
  
  patient_result <- expand.grid(simulation = i, method = methods,
                                patient_id = seq_len(nrow(fixed_patients)),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  patient_result$p_true <- p_true_pat[patient_result$patient_id]
  patient_result$p_estimated <- NA_real_
  
  # Store validation performance and all fixed-patient predictions together
  save_method <- function(label, p_val, p_patient, varsel, option = NA_real_) {
    idx <- match(label, methods)
    scores <- as.numeric(measures(yval, as.numeric(p_val), p_true_val))
    method_result[idx, ] <<- c(prev, auc, ndev, label, scores, as.integer(varsel), option)
    patient_result$p_estimated[patient_result$method == label] <<- as.numeric(p_patient)
  }
  
  # Full maximum likelihood logistic regression
  fit <- glm(y ~ ., data = data.dev, family = "binomial")
  p_val <- as.numeric(predict(fit, newdata = data.val, type = "response"))
  p_pat <- as.numeric(predict(fit, newdata = fixed_x, type = "response"))
  save_method("MLE", p_val, p_pat, rep(1L, n.para))
  
  # Univariable selection at 5%
  unisum_05 <- unilogit(x, y, 0.05)
  varsel_uni_05 <- unisum_05$varsel_uni
  unimodel_05 <- unisum_05$uni
  selected <- which(as.integer(varsel_uni_05) == 1L)
  Xval_uni <- cbind(1, as.matrix(xval[, selected, drop = FALSE]))
  Xpat_uni <- cbind(1, as.matrix(fixed_x[, selected, drop = FALSE]))
  p_val <- plogis(as.vector(Xval_uni %*% coef(unimodel_05)))
  p_pat <- plogis(as.vector(Xpat_uni %*% coef(unimodel_05)))
  save_method("UVS-5%", p_val, p_pat, varsel_uni_05, 0.05)
  
  # Univariable selection at 15%
  unisum_15 <- unilogit(x, y, 0.15)
  varsel_uni_15 <- unisum_15$varsel_uni
  unimodel_15 <- unisum_15$uni
  selected <- which(as.integer(varsel_uni_15) == 1L)
  Xval_uni <- cbind(1, as.matrix(xval[, selected, drop = FALSE]))
  Xpat_uni <- cbind(1, as.matrix(fixed_x[, selected, drop = FALSE]))
  p_val <- plogis(as.vector(Xval_uni %*% coef(unimodel_15)))
  p_pat <- plogis(as.vector(Xpat_uni %*% coef(unimodel_15)))
  save_method("UVS-15%", p_val, p_pat, varsel_uni_15, 0.15)
  
  # Backward elimination at 5%
  back_05 <- backward_pvalue(data.dev, "y", 0.05)
  varsel_back_05 <- back_05$varsel_back
  backmodel_05 <- back_05$backmodel
  p_val <- as.numeric(predict(backmodel_05, newdata = data.val, type = "response"))
  p_pat <- as.numeric(predict(backmodel_05, newdata = fixed_x, type = "response"))
  save_method("BE-5%", p_val, p_pat, varsel_back_05, 0.05)
  
  # Backward elimination at 15%
  back_15 <- backward_pvalue(data.dev, "y", 0.15)
  varsel_back_15 <- back_15$varsel_back
  backmodel_15 <- back_15$backmodel
  p_val <- as.numeric(predict(backmodel_15, newdata = data.val, type = "response"))
  p_pat <- as.numeric(predict(backmodel_15, newdata = fixed_x, type = "response"))
  save_method("BE-15%", p_val, p_pat, varsel_back_15, 0.15)
  
  # Standard LASSO at lambda.min
  lasso <- glmnet::cv.glmnet(x_mat, y, alpha = 1, family = "binomial", type.measure = "deviance")
  lambda_min <- lasso$lambda.min
  p_val <- as.numeric(predict(lasso, newx = xval_mat, s = lambda_min, type = "response"))
  p_pat <- as.numeric(predict(lasso, newx = fixed_mat, s = lambda_min, type = "response"))
  varsel_min <- as.integer(as.numeric(coef(lasso, s = lambda_min))[-1] != 0)
  save_method("LS-min", p_val, p_pat, varsel_min, lambda_min)
  
  # Modified LASSO using bootstrap-averaged coefficients
  nfolds <- 10
  f <- nfolds / (nfolds - 1) - 1
  mod_lasso <- mod_penal_ave_foreach(x = x_mat, y = y, bn = 20, method = "lasso",
                                     f = f, parallel = FALSE, nfolds = nfolds, boot = TRUE)
  beta_mod <- as.numeric(mod_lasso$beta.boot)
  if (length(beta_mod) != n.para + 1L) stop("Invalid modified LASSO coefficients")
  p_val <- plogis(as.vector(cbind(1, xval_mat) %*% beta_mod))
  p_pat <- plogis(as.vector(cbind(1, fixed_mat) %*% beta_mod))
  varsel_mod <- as.integer(beta_mod[-1] != 0)
  save_method("LS-mod", p_val, p_pat, varsel_mod, as.numeric(mod_lasso$lambda.boot)[1])
  
  # Adaptive LASSO with Ridge-based penalty weights
  ridge_init <- glmnet::cv.glmnet(x_mat, y, alpha = 0, family = "binomial",
                                  type.measure = "deviance", nfolds = 10)
  beta_initial <- as.numeric(coef(ridge_init, s = "lambda.min"))[-1]
  adaptive_weights <- 1 / (abs(beta_initial) + 1e-6)
  adaptive_lasso <- glmnet::cv.glmnet(x_mat, y, alpha = 1, family = "binomial",
                                      penalty.factor = adaptive_weights,
                                      type.measure = "deviance", nfolds = 10)
  p_val <- as.numeric(predict(adaptive_lasso, newx = xval_mat,
                              s = "lambda.min", type = "response"))
  p_pat <- as.numeric(predict(adaptive_lasso, newx = fixed_mat,
                              s = "lambda.min", type = "response"))
  varsel_adaptive <- as.integer(as.numeric(coef(adaptive_lasso, s = "lambda.min"))[-1] != 0)
  save_method("ALASSO", p_val, p_pat, varsel_adaptive, adaptive_lasso$lambda.min)
  
  # Harrell step-down approximation to MLE linear predictor
  eta_mle <- as.numeric(predict(fit, type = "link"))
  sd_mle <- prune_predictor_harrell(x, eta_mle, r2_threshold = 0.99)
  p_val <- predict_harrell_stepdown(sd_mle, xval)
  p_pat <- predict_harrell_stepdown(sd_mle, fixed_x)
  save_method("SD-MLE", p_val, p_pat, sd_mle$kept, sd_mle$r2)
  
  # Harrell step-down approximation to LASSO linear predictor
  eta_lasso <- as.numeric(predict(lasso, newx = x_mat, s = lambda_min, type = "link"))
  sd_lasso <- prune_predictor_harrell(x[,varsel_min == 1], eta_lasso, r2_threshold = 0.99)
  p_val <- predict_harrell_stepdown(sd_lasso, xval)
  p_pat <- predict_harrell_stepdown(sd_lasso, fixed_x)
  save_method("SD-LS-min", p_val, p_pat, sd_lasso$kept, sd_lasso$r2)
  
  # Harrell step-down approximation to modified LASSO linear predictor
  eta_mod <- as.vector(cbind(1, x_mat) %*% beta_mod)
  sd_mod <- prune_predictor_harrell(x[,varsel_mod == 1], eta_mod, r2_threshold = 0.99)
  p_val <- predict_harrell_stepdown(sd_mod, xval)
  p_pat <- predict_harrell_stepdown(sd_mod, fixed_x)
  save_method("SD-LS-mod", p_val, p_pat, sd_mod$kept, sd_mod$r2)
  
  # Return performance and fixed-patient predictions
  colnames(method_result) <- c("prevalence", "anticipated c-stat", "ndev", "method",
                               "calibration slope", "calibration in the large", "auc",
                               "brier", "mape", paste0("varsel", seq_len(n.para)), "option")
  method_result <- as.data.frame(method_result, stringsAsFactors = FALSE)
  numeric_cols <- setdiff(names(method_result), "method")
  method_result[numeric_cols] <- lapply(method_result[numeric_cols], as.numeric)
  method_result$simulation <- i
  
  rownames(patient_result) <- NULL
  patient_result$simulation <- as.integer(patient_result$simulation)
  patient_result$patient_id <- as.integer(patient_result$patient_id)
  patient_result$p_true <- as.numeric(patient_result$p_true)
  patient_result$p_estimated <- as.numeric(patient_result$p_estimated)
  
  list(performance = method_result, patient_predictions = patient_result)
}
