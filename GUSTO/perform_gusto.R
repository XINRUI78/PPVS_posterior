library(doParallel)
library(foreach)

perform_gusto <- function(order, data, n.para, n.restrict = NULL,
                    thresholds = seq(0.01, 0.50, by = 0.01)) {
  
  ndev <- length(order)
  
  # Development dataset
  data.dev <- data[order, ]
  x <- data.dev[, -1]
  y <- data.dev[, 1]
  
  # Validation dataset
  data.val <- data[-order, ]
  
  yval <- data.val[, 1]
  
  ##########################################################
  # Initialise storage
  ##########################################################
  
  n_methods <- if (is.null(n.restrict)) 14 else 17
  
  # ndev + method + 4 performance measures
  # + n.para variable-selection indicators + option
  method_result <- matrix(NA, nrow = n_methods, ncol = 7 + n.para)
  
  # Store net benefit separately
  nb_list <- list()
  
  
  ##########################################################
  # Method 0: Full MLE
  ##########################################################
  
  fit <- glm(y ~ ., data = data.dev, family = "binomial")
  xval <- model.matrix(fit, data = data.val)[,-1]
  eta_val <- as.matrix(cbind(1, xval)) %*% coef(fit)
  p_val <- as.vector(1 / (1 + exp(-eta_val)))
  
  method_result[1, ] <- c(
    ndev, 0,
    measures(yval, p_val),
    rep(1, n.para),
    NA
  )
  
  nb_list[["MLE"]] <- calc_net_benefit(yval, p_val, thresholds, "MLE")
  
  
  ##########################################################
  # Method 1: Backward elimination p < 0.05
  ##########################################################
  
  p_threshold <- 0.05
  back_05 <- backward_pvalue(data.dev, "y", p_threshold)
  varsel_back_05 <- back_05$varsel_back
  backmodel_05 <- back_05$backmodel
  
  back_eta_05 <- as.matrix(
    cbind(1, xval[, varsel_back_05 == 1, drop = FALSE])
  ) %*% coef(backmodel_05)
  
  back_p_05 <- as.vector(1 / (1 + exp(-back_eta_05)))
  
  method_result[2, ] <- c(
    ndev, 1,
    measures(yval, back_p_05),
    varsel_back_05,
    p_threshold
  )
  
  nb_list[["BE_05"]] <- calc_net_benefit(yval, back_p_05, thresholds, "BE_05")
  
  
  ##########################################################
  # Method 2: Backward elimination p < 0.15
  ##########################################################
  
  p_threshold <- 0.15
  back_15 <- backward_pvalue(data.dev, "y", p_threshold)
  varsel_back_15 <- back_15$varsel_back
  backmodel_15 <- back_15$backmodel
  
  back_eta_15 <- as.matrix(
    cbind(1, xval[, varsel_back_15 == 1, drop = FALSE])
  ) %*% coef(backmodel_15)
  
  back_p_15 <- as.vector(1 / (1 + exp(-back_eta_15)))
  
  method_result[3, ] <- c(
    ndev, 2,
    measures(yval, back_p_15),
    varsel_back_15,
    p_threshold
  )
  
  nb_list[["BE_15"]] <- calc_net_benefit(yval, back_p_15, thresholds, "BE_15")
  
  
  ##########################################################
  # Method 3: Univariable selection p < 0.05
  ##########################################################
  
  p_threshold <- 0.05
  unisum_05 <- unilogit(x, y, p_threshold)
  varsel_uni_05 <- unisum_05$varsel_uni
  unimodel_05 <- unisum_05$uni
  
  uni_eta_05 <- as.matrix(
    cbind(1, xval[, varsel_uni_05 == 1, drop = FALSE])
  ) %*% coef(unimodel_05)
  
  uni_p_05 <- as.vector(1 / (1 + exp(-uni_eta_05)))
  
  method_result[4, ] <- c(
    ndev, 3,
    measures(yval, uni_p_05),
    varsel_uni_05,
    p_threshold
  )
  
  nb_list[["UVS_05"]] <- calc_net_benefit(yval, uni_p_05, thresholds, "UVS_05")
  
  
  ##########################################################
  # Method 4: Univariable selection p < 0.15
  ##########################################################
  
  p_threshold <- 0.15
  unisum_15 <- unilogit(x, y, p_threshold)
  varsel_uni_15 <- unisum_15$varsel_uni
  unimodel_15 <- unisum_15$uni
  
  uni_eta_15 <- as.matrix(
    cbind(1, xval[, varsel_uni_15 == 1, drop = FALSE])
  ) %*% coef(unimodel_15)
  
  uni_p_15 <- as.vector(1 / (1 + exp(-uni_eta_15)))
  
  method_result[5, ] <- c(
    ndev, 4,
    measures(yval, uni_p_15),
    varsel_uni_15,
    p_threshold
  )
  
  nb_list[["UVS_15"]] <- calc_net_benefit(yval, uni_p_15, thresholds, "UVS_15")
  
  
  ##########################################################
  # Method 5: UVS 5% + BE 5%
  ##########################################################
  
  p_threshold <- 0.05
  data.dev2 <- data.dev[, c(1, 1 + which(varsel_uni_05 == 1)), drop = FALSE]
  uni05_back <- backward_pvalue(data.dev2, "y", p_threshold)
  
  varsel_uni05_back <- varsel_uni_05
  varsel_uni05_back[varsel_uni05_back == 1] <- uni05_back$varsel_back
  uni05_back_model <- uni05_back$backmodel
  
  uni05_back_eta <- as.matrix(
    cbind(1, xval[, varsel_uni05_back == 1, drop = FALSE])
  ) %*% coef(uni05_back_model)
  
  uni05_back_p <- as.vector(1 / (1 + exp(-uni05_back_eta)))
  
  method_result[6, ] <- c(
    ndev, 5,
    measures(yval, uni05_back_p),
    varsel_uni05_back,
    0.05
  )
  
  nb_list[["UVS05_BE05"]] <- calc_net_benefit(yval, uni05_back_p, thresholds, "UVS05_BE05")
  
  
  ##########################################################
  # Method 6: UVS 15% + BE 5%
  ##########################################################
  
  p_threshold <- 0.05
  data.dev2 <- data.dev[, c(1, 1 + which(varsel_uni_15 == 1)), drop = FALSE]
  uni15_back <- backward_pvalue(data.dev2, "y", p_threshold)
  
  varsel_uni15_back <- varsel_uni_15
  varsel_uni15_back[varsel_uni15_back == 1] <- uni15_back$varsel_back
  uni15_back_model <- uni15_back$backmodel
  
  uni15_back_eta <- as.matrix(
    cbind(1, xval[, varsel_uni15_back == 1, drop = FALSE])
  ) %*% coef(uni15_back_model)
  
  uni15_back_p <- as.vector(1 / (1 + exp(-uni15_back_eta)))
  
  method_result[7, ] <- c(
    ndev, 6,
    measures(yval, uni15_back_p),
    varsel_uni15_back,
    0.15
  )
  
  nb_list[["UVS15_BE05"]] <- calc_net_benefit(yval, uni15_back_p, thresholds, "UVS15_BE05")
  
  
  ##########################################################
  # Methods 7 and 8: LASSO lambda.min and lambda.1se
  ##########################################################
  
  lasso <- glmnet::cv.glmnet(
    as.matrix(x), y, alpha = 1,
    family = "binomial", type.measure = "deviance"
  )
  
  lambda_min <- lasso$lambda.min
  lambda_1se <- lasso$lambda.1se
  
  lassomin_p <- as.vector(
    predict(lasso, as.matrix(xval), s = lambda_min, type = "response")
  )
  
  lasso1se_p <- as.vector(
    predict(lasso, as.matrix(xval), s = lambda_1se, type = "response")
  )
  
  varsel_min <- ifelse(
    as.numeric(coef(lasso, s = lambda_min)[-1]) != 0, 1, 0
  )
  
  varsel_1se <- ifelse(
    as.numeric(coef(lasso, s = lambda_1se)[-1]) != 0, 1, 0
  )
  
  method_result[8, ] <- c(
    ndev, 7,
    measures(yval, lassomin_p),
    varsel_min,
    lambda_min
  )
  
  method_result[9, ] <- c(
    ndev, 8,
    measures(yval, lasso1se_p),
    varsel_1se,
    lambda_1se
  )
  
  nb_list[["LS_min"]] <- calc_net_benefit(yval, lassomin_p, thresholds, "LS_min")
  
  nb_list[["LS_1se"]] <- calc_net_benefit(yval, lasso1se_p, thresholds, "LS_1se")
  
  
  ##########################################################
  # Method 9: LASSO lambda.min + MLE
  ##########################################################
  
  data.s2 <- data.dev[, c(1, 1 + which(varsel_min == 1)), drop = FALSE]
  fit <- glm(y ~ ., data = data.s2, family = "binomial")
  
  if (sum(varsel_min) == 0) {
    eta_min_mle <- rep(coef(fit)[1], nrow(xval))
  } else {
    eta_min_mle <- as.matrix(
      cbind(1, xval[, varsel_min == 1, drop = FALSE])
    ) %*% coef(fit)
  }
  
  p_min_mle <- as.vector(1 / (1 + exp(-eta_min_mle)))
  
  method_result[10, ] <- c(
    ndev, 9,
    measures(yval, p_min_mle),
    varsel_min,
    lambda_min
  )
  
  nb_list[["LSmin_MLE"]] <- calc_net_benefit(yval, p_min_mle, thresholds, "LSmin_MLE")
  
  
  ##########################################################
  # Method 10: LASSO lambda.1se + MLE
  ##########################################################
  
  data.s2 <- as.data.frame(
    data.dev[, c(1, 1 + which(varsel_1se == 1)), drop = FALSE]
  )
  
  fit <- glm(y ~ ., data = data.s2, family = binomial())
  
  if (sum(varsel_1se) == 0) {
    eta_1se_mle <- rep(coef(fit)[1], nrow(xval))
  } else {
    eta_1se_mle <- as.matrix(
      cbind(1, xval[, varsel_1se == 1, drop = FALSE])
    ) %*% coef(fit)
  }
  
  p_1se_mle <- as.vector(1 / (1 + exp(-eta_1se_mle)))
  
  method_result[11, ] <- c(
    ndev, 10,
    measures(yval, p_1se_mle),
    varsel_1se,
    lambda_1se
  )
  
  nb_list[["LS1se_MLE"]] <- calc_net_benefit(yval, p_1se_mle, thresholds, "LS1se_MLE")
  
  
  ##########################################################
  # Method 11: LASSO lambda.min + BE
  ##########################################################
  
  p_threshold <- 0.05
  data.s2 <- data.dev[, c(1, 1 + which(varsel_min == 1)), drop = FALSE]
  
  if (sum(varsel_min) == 0) {
    fit <- glm(y ~ 1, data = data.dev, family = binomial())
    lasso_min_back_eta <- rep(coef(fit)[1], nrow(xval))
    varsel_min_back <- varsel_min
  } else {
    lassomin_back <- backward_pvalue(data.s2, "y", p_threshold)
    varsel_min_back <- varsel_min
    varsel_min_back[varsel_min_back == 1] <- lassomin_back$varsel_back
    lasso_min_back_model <- lassomin_back$backmodel
    
    if (sum(varsel_min_back) == 0) {
      lasso_min_back_eta <- rep(coef(lasso_min_back_model)[1], nrow(xval))
    } else {
      lasso_min_back_eta <- as.matrix(
        cbind(1, xval[, varsel_min_back == 1, drop = FALSE])
      ) %*% coef(lasso_min_back_model)
    }
  }
  
  lasso_min_back_p <- as.vector(1 / (1 + exp(-lasso_min_back_eta)))
  
  method_result[12, ] <- c(
    ndev, 11,
    measures(yval, lasso_min_back_p),
    varsel_min_back,
    lambda_min
  )
  
  nb_list[["LSmin_BE05"]] <- calc_net_benefit(yval, lasso_min_back_p, thresholds, "LSmin_BE05")
  
  
  ##########################################################
  # Method 12: LASSO lambda.1se + BE
  ##########################################################
  
  p_threshold <- 0.05
  data.s2 <- as.data.frame(
    data.dev[, c(1, 1 + which(varsel_1se == 1)), drop = FALSE]
  )
  
  if (sum(varsel_1se) == 0) {
    fit <- glm(y ~ 1, data = data.dev, family = binomial())
    lasso1se_back_eta <- rep(coef(fit)[1], nrow(xval))
    varsel_lasso1se_back <- varsel_1se
  } else {
    lasso_1se_back <- backward_pvalue(data.s2, "y", p_threshold)
    varsel_lasso1se_back <- varsel_1se
    varsel_lasso1se_back[varsel_lasso1se_back == 1] <- lasso_1se_back$varsel_back
    lasso1se_back_model <- lasso_1se_back$backmodel
    
    if (sum(varsel_lasso1se_back) == 0) {
      lasso1se_back_eta <- rep(coef(lasso1se_back_model)[1], nrow(xval))
    } else {
      lasso1se_back_eta <- as.matrix(
        cbind(1, xval[, varsel_lasso1se_back == 1, drop = FALSE])
      ) %*% coef(lasso1se_back_model)
    }
  }
  
  lasso1se_back_p <- as.vector(1 / (1 + exp(-lasso1se_back_eta)))
  
  method_result[13, ] <- c(
    ndev, 12,
    measures(yval, lasso1se_back_p),
    varsel_lasso1se_back,
    lambda_1se
  )
  
  nb_list[["LS1se_BE05"]] <- calc_net_benefit(yval, lasso1se_back_p, thresholds, "LS1se_BE05")
  
  
  ##########################################################
  # Method 13: Modified LASSO
  ##########################################################
  
  nfolds <- 10
  f <- nfolds / (nfolds - 1) - 1
  
  mod_lasso <- mod_penal_ave_foreach(
    x = as.matrix(x), y = y, bn = 20,
    method = "lasso", f = f,
    parallel = FALSE, nfolds = nfolds, boot = TRUE
  )
  
  eta_mod_lasso <- as.matrix(cbind(1, xval)) %*% mod_lasso$beta.boot
  p_mod_lasso <- as.vector(1 / (1 + exp(-eta_mod_lasso)))
  
  varsel_mod_lasso <- ifelse(
    as.numeric(mod_lasso$beta.boot)[-1] != 0, 1, 0
  )
  
  method_result[14, ] <- c(
    ndev, 13,
    measures(yval, p_mod_lasso),
    varsel_mod_lasso,
    as.numeric(mod_lasso$lambda.boot)
  )
  
  nb_list[["LS_mod"]] <- calc_net_benefit(yval, p_mod_lasso, thresholds, "LS_mod")
  
  
  ##########################################################
  # Restricted methods
  ##########################################################
  
  if (!is.null(n.restrict)) {
    
    # Method 14: Univariable ranking
    unirank_result <- unirank(x, y, n.restrict)
    varsel_uni <- unirank_result$varsel_uni
    unimodel <- unirank_result$model
    
    uni_eta <- as.matrix(
      cbind(1, xval[, varsel_uni == 1, drop = FALSE])
    ) %*% coef(unimodel)
    
    uni_p <- as.vector(1 / (1 + exp(-uni_eta)))
    
    method_result[15, ] <- c(
      ndev, 14,
      measures(yval, uni_p),
      varsel_uni,
      NA
    )
    
    nb_list[["UVS_restricted"]] <- calc_net_benefit(yval, uni_p, thresholds, "UVS_restricted")
    
    
    # Method 15: Backward elimination ranking
    berank_result <- berank(x, y, n.restrict)
    varsel_be <- berank_result$varsel_be
    bemodel <- berank_result$model
    
    be_eta <- as.matrix(
      cbind(1, xval[, varsel_be == 1, drop = FALSE])
    ) %*% coef(bemodel)
    
    be_p <- as.vector(1 / (1 + exp(-be_eta)))
    
    method_result[16, ] <- c(
      ndev, 15,
      measures(yval, be_p),
      varsel_be,
      NA
    )
    
    nb_list[["BE_restricted"]] <- calc_net_benefit(yval, be_p, thresholds, "BE_restricted")
    
    
    # Method 16: Restricted LASSO
    lasso_exact_result <- lasso_exact(
      x, y, xval, n.restrict,
      max_attempts = 10, initial_nlambda = 100
    )
    
    varsel_lasso <- lasso_exact_result$varsel_lasso
    lambda <- lasso_exact_result$lambda
    lasso_p <- lasso_exact_result$lasso_p
    
    method_result[17, ] <- c(
      ndev, 16,
      measures(yval, lasso_p),
      varsel_lasso,
      lambda
    )
    
    nb_list[["LASSO_restricted"]] <- calc_net_benefit(yval, lasso_p, thresholds, "LASSO_restricted")
  }
  
  
  ##########################################################
  # Treat-all and treat-none
  ##########################################################
  
  prevalence_val <- mean(yval)
  
  nb_list[["Treat_all"]] <- data.frame(
    method = "Treat_all",
    threshold = thresholds,
    net_benefit = prevalence_val -
      (1 - prevalence_val) * thresholds / (1 - thresholds)
  )
  
  nb_list[["Treat_none"]] <- data.frame(
    method = "Treat_none",
    threshold = thresholds,
    net_benefit = 0
  )
  
  
  ##########################################################
  # Column names
  ##########################################################
  
  colnames(method_result) <- c(
    "ndev",
    "method",
    "calibration slope",
    "calibration in the large",
    "auc",
    "Brier score",
    paste0("varsel", 1:n.para),
    "option"
  )
  
  
  ##########################################################
  # Net benefit long format
  ##########################################################
  
  net_benefit_long <- do.call(rbind, nb_list)
  rownames(net_benefit_long) <- NULL
  
  
  ##########################################################
  # Return both
  ##########################################################
  
  return(list(
    method_result = method_result,
    net_benefit = net_benefit_long
  ))
}
