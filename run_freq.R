
# ============================================================
# run_freq.R
# Nine frequentist methods
# One replication per call
# ============================================================

library(glmnet)

run_freq <- function(
    i, scenario_id, ndev, n.para,
    beta0, beta, fixed_patients,
    data_dir = "simulation_data"
) {

  # ----------------------------------------------------------
  # 1. Read the existing simulation dataset
  # ----------------------------------------------------------

  data_file <- simulation_file(
    scenario_id, ndev, i, data_dir
  )

  sim_data <- readRDS(data_file)
  data.dev <- sim_data$data.dev
  data.val <- sim_data$data.val

  predictor_names <- paste0("X", seq_len(n.para))
  x <- data.dev[, predictor_names, drop = FALSE]
  y <- as.numeric(data.dev$y)

  xval <- data.val[, predictor_names, drop = FALSE]
  yval <- as.numeric(data.val$y)

  fixed_x <- fixed_patients[
    , predictor_names, drop = FALSE
  ]

  patient_names <- c(
    "low1", "low2",
    "medium1", "medium2",
    "high1", "high2"
  )

  # True probabilities for validation observations
  p_true_val <- plogis(as.numeric(beta0 + as.matrix(xval) %*% as.numeric(beta)))

  # True probabilities for the fixed patients
  p_true_pat <- as.numeric(fixed_patients$p_true)

  # Keep the original performance matrix structure (9 methods).
  prev <- 0.3
  auc <- 0.8

  method_result <- matrix(NA_character_, nrow = 9, ncol = 10 + n.para)
  # ========================================================
    # Initialise six-patient probability results
    # ========================================================

    patient_id <- seq_len(nrow(fixed_patients))
    fixed_p_true <- p_true_pat

    patient_result <- expand.grid(
      simulation = i,
      method = c(
        "MLE", "UVS-5%", "UVS-15%",
        "BE-5%", "BE-15%",
        "LS-min", "step-min", "LS-mod", "step-mod"
      ),
      patient_id = patient_id,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )


    patient_result$p_true <- fixed_p_true[patient_result$patient_id]

    patient_result$p_estimated <- NA_real_

  # ----------------------------------------------------------
  # Method 1: Full MLE
  # ----------------------------------------------------------
   
    fit <- glm(y ~ ., data = data.dev, family = 'binomial')
    eta_val <- as.matrix(cbind(1,xval))%*%coef(fit)
    p_val <- as.vector(1/(1+exp(-eta_val)))
    method_result[1,] <- c(prev, auc, ndev, "MLE", measures(yval, p_val, p_true_val), rep(1, n.para), NA)

   # Six fixed-patient predictions
      p_patient_mle <- as.numeric(
        predict(
          fit,
          newdata = fixed_x,
          type = "response"
        )
      )

      patient_result$p_estimated[
        patient_result$method == "MLE"
      ] <- p_patient_mle
    
  # ----------------------------------------------------------
  # Methods 2 and 3: Univariable selection
  # ----------------------------------------------------------
  p_threshold <- 0.05
  unisum_05 <- unilogit(x, y, p_threshold)
  varsel_uni_05 <- unisum_05$varsel_uni
  unimodel_05 <- unisum_05$uni
  uni_eta_05 <- as.matrix(cbind(1, xval[, varsel_uni_05 == 1])) %*% coef(unimodel_05)
  uni_p_05 <- as.vector(1 / (1 + exp(-uni_eta_05)))
  method_result[2,] <- c(prev, auc, ndev, "UVS-5%", measures(yval, uni_p_05, p_true_val), varsel_uni_05, p_threshold)

  p_patient_uni5 <- as.numeric(
        predict(unimodel_05, newdata = fixed_x, type = "response")
      )

      patient_result$p_estimated[
        patient_result$method == "UVS-5%"
      ] <- p_patient_uni5
  
  p_threshold <- 0.15
  unisum_15 <- unilogit(x, y, p_threshold)
  varsel_uni_15 <- unisum_15$varsel_uni
  unimodel_15 <- unisum_15$uni
  uni_eta_15 <- as.matrix(cbind(1, xval[, varsel_uni_15 == 1])) %*% coef(unimodel_15)
  uni_p_15 <- as.vector(1 / (1 + exp(-uni_eta_15)))
  method_result[3, ] <- c(prev, auc, ndev, "UVS-15%", measures(yval, uni_p_15, p_true_val), varsel_uni_15, p_threshold)

  p_patient_uni15 <- as.numeric(
        predict(unimodel_15, newdata = fixed_x, type = "response")
      )

      patient_result$p_estimated[
        patient_result$method == "UVS-15%"
      ] <- p_patient_uni15   

  # ----------------------------------------------------------
  # Methods 4 and 5: Backward elimination
  # ----------------------------------------------------------
  p_threshold = 0.05
  back_05 <- backward_pvalue(data.dev, "y", p_threshold) 
  varsel_back_05 <- back_05$varsel_back
  backmodel_05 <- back_05$backmodel
  back_eta_05 <- as.matrix(cbind(1,xval[,varsel_back_05 == 1]))%*%coef(backmodel_05)
  back_p_05 <- as.vector(1/(1+exp(-back_eta_05)))
  method_result[4,] <- c(prev, auc, ndev, "BE-5%", measures(yval, back_p_05, p_true_val), varsel_back_05, p_threshold)
 
  p_patient_be5 <- as.numeric(
        predict(backmodel_05, newdata = fixed_x, type = "response")
      )

      patient_result$p_estimated[
        patient_result$method == "BE-5%"
      ] <- p_patient_be5   

  p_threshold = 0.15
  back_15 <- backward_pvalue(data.dev, "y", p_threshold) 
  varsel_back_15 <- back_15$varsel_back
  backmodel_15 <- back_15$backmodel
  back_eta_15 <- as.matrix(cbind(1,xval[,varsel_back_15 == 1]))%*%coef(backmodel_15)
  back_p_15 <- as.vector(1/(1+exp(-back_eta_15)))
  method_result[5,] <- c(prev, auc, ndev, "BE-15%", measures(yval, back_p_15, p_true_val), varsel_back_15, p_threshold)

  p_patient_be15 <- as.numeric(
        predict(backmodel_15, newdata = fixed_x, type = "response")
      )

      patient_result$p_estimated[
        patient_result$method == "BE-15%"
      ] <- p_patient_be15  

  # ----------------------------------------------------------
  # Method 6: Standard LASSO using lambda.min
  # ----------------------------------------------------------

  lasso <- glmnet::cv.glmnet(as.matrix(x), y, alpha = 1, family = "binomial", type.measure = "deviance")
  lambda_min <- lasso$lambda.min
  # validate the fitted models on the validation dataset
  lassomin_p <- as.vector(predict(lasso, as.matrix(xval), s = lambda_min, type="response"))
  varsel_min <- ifelse(as.numeric(coef(lasso, s = lambda_min)[-1]) != 0, 1, 0)
  method_result[6,] <- c(prev, auc, ndev, "LS-min", measures(yval, lassomin_p, p_true_val), varsel_min, lambda_min)

  p_patient_lasso <- as.numeric(
        predict(
          lasso,
          newx = as.matrix(fixed_x),
          s = lambda_min,
          type = "response"
        )
      )

      patient_result$p_estimated[
        patient_result$method == "LS-min"
      ] <- p_patient_lasso

  
# ============================================================
# Method 7: LASSO min with R-squared pruning
# ============================================================

  eta_lasso <- as.numeric(
    predict(lasso, newx = as.matrix(x), s = "lambda.min", type = "link")
  )

  pruned <- prune_lasso_predictor(
    x = x,
    eta_lasso = eta_lasso,
    r2_threshold = 0.99
  )

  beta_pruned <- pruned$coefficients
  varsel_stepmin <- pruned$kept

  p_val <- plogis(as.numeric(model.matrix(~ ., data = xval)[, names(beta_pruned), drop = FALSE] %*% beta_pruned))
  method_result[7,] <- c(prev, auc, ndev, "step-min", measures(yval, p_val, p_true_val), varsel_stepmin, NA)

  p_patient_stepmin <- plogis(as.numeric(model.matrix(~ ., data = fixed_x)[, names(beta_pruned), drop = FALSE] %*% beta_pruned))
  
  patient_result$p_estimated[
        patient_result$method == "step-min"
      ] <- p_patient_stepmin


  # ----------------------------------------------------------
  # Method 7: Modified LASSO
  # ----------------------------------------------------------

  nfolds <- 10
  f <- nfolds / (nfolds - 1) - 1
  mod_lasso <- mod_penal_ave_foreach(x = as.matrix(x), y = y, bn = 20, method = "lasso", f = f,
                  parallel = FALSE, nfolds = nfolds, boot = TRUE)
  eta_mod_lasso <- as.matrix(cbind(1,xval))%*% mod_lasso$beta.boot
  p_mod_lasso <- as.vector(1/(1+exp(-eta_mod_lasso)))
  varsel_mod_lasso <- ifelse(as.numeric(mod_lasso$beta.boot)[-1] != 0, 1, 0)
  method_result[8, ] <- c(prev, auc, ndev, "LS-mod", measures(yval, p_mod_lasso, p_true_val), varsel_mod_lasso, as.numeric(mod_lasso$lambda.boot))
 
  eta_fix <- as.matrix(cbind(1,fixed_x))%*% mod_lasso$beta.boot
  p_fix <- as.vector(1/(1+exp(-eta_fix)))

      patient_result$p_estimated[
        patient_result$method == "LS-mod"
      ] <- p_fix
      #########
    eta_lasso <- as.numeric(
    as.matrix(cbind(1, x)) %*% mod_lasso$beta.boot
  )

  pruned <- prune_lasso_predictor(
    x = x,
    eta_lasso = eta_lasso,
    r2_threshold = 0.99
  )

  beta_pruned <- pruned$coefficients
  varsel_stepmod <- pruned$kept

  p_val <- plogis(as.numeric(model.matrix(~ ., data = xval)[, names(beta_pruned), drop = FALSE] %*% beta_pruned))
  method_result[9,] <- c(prev, auc, ndev, "step-mod", measures(yval, p_val, p_true_val), varsel_stepmod, NA)

  p_patient_stepmod <- plogis(as.numeric(model.matrix(~ ., data = fixed_x)[, names(beta_pruned), drop = FALSE] %*% beta_pruned))
  
  patient_result$p_estimated[
        patient_result$method == "step-mod"
      ] <- p_patient_stepmod

  # ----------------------------------------------------------
  # Return both tasks
  # ----------------------------------------------------------
  patient_predictions <- patient_result
  rownames(patient_predictions) <- NULL

  # Ensure numerical patient columns
  patient_predictions$simulation <- as.integer(patient_predictions$simulation)
  patient_predictions$patient_id <- as.integer(patient_predictions$patient_id)
  patient_predictions$observed_y <- as.numeric(patient_predictions$observed_y)
  patient_predictions$p_true <- as.numeric(patient_predictions$p_true)
  patient_predictions$p_estimated <- as.numeric(patient_predictions$p_estimated)

  colnames(method_result) <- c("prevalence",
                         "anticipated c-stat",
                         "ndev",
                         "method",
                         "calibration slope",
                         "calibration in the large",
                         "auc",
                         "brier",
                         "mape",
                         paste0("varsel", 1:n.para),
                         "option")

  method_result <- as.data.frame(method_result, stringsAsFactors = FALSE)
  numeric_cols <- setdiff(names(method_result), "method")
  method_result[numeric_cols] <- lapply(method_result[numeric_cols], as.numeric)
  method_result$simulation <- i


  # ==========================================================
  # Return both outputs
  # ==========================================================

  list(
    performance = method_result,
    patient_predictions = patient_predictions
  )
}

