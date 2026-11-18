library(doParallel)
library(foreach)

perform_gusto <- function(order, data, n.para, n.restrict = NULL,
                          thresholds = seq(0.01, 0.50, by = 0.01)) {
  ndev <- length(order)
  data.dev <- data[order, , drop = FALSE]
  data.val <- data[-order, , drop = FALSE]
  x <- data.dev[, -1, drop = FALSE]
  y <- data.dev[, 1]
  xval <- data.val[, -1, drop = FALSE]
  yval <- data.val[, 1]
  predictors <- names(x)
  if (length(predictors) != n.para) stop("n.para does not match predictors")

  method_result <- matrix(NA_real_, nrow = 13, ncol = 7 + n.para)
  nb_list <- list()
  store <- function(row, p, varsel, option = NA_real_, name) {
    method_result[row, ] <<- c(ndev, row, measures(yval, p), varsel, option)
    nb_list[[name]] <<- calc_net_benefit(yval, p, thresholds, name)
  }
  predict_glm <- function(model) as.vector(predict(model, newdata = data.val,
                                                   type = "response"))

  fit_full <- glm(reformulate(predictors, response = names(data.dev)[1]),
                  data = data.dev, family = binomial())
  eta_mle <- as.vector(cbind(1, as.matrix(x)) %*% coef(fit_full))
  p_val <- plogis(as.vector(cbind(1, as.matrix(xval)) %*% coef(fit_full)))
  store(1, p_val, rep(1, n.para), name = "MLE")

  back_05 <- backward_pvalue(data.dev, "y", 0.05)
  back_p_05 <- predict_glm(back_05$backmodel)
  store(2, back_p_05, back_05$varsel_back, 0.05, "BE_05")

  back_15 <- backward_pvalue(data.dev, "y", 0.15)
  back_p_15 <- predict_glm(back_15$backmodel)
  store(3, back_p_15, back_15$varsel_back, 0.15, "BE_15")

  unisum_05 <- unilogit(x, y, 0.05)
  uni_p_05 <- predict_glm(unisum_05$uni)
  store(4, uni_p_05, unisum_05$varsel_uni, 0.05, "UVS_05")

  unisum_15 <- unilogit(x, y, 0.15)
  uni_p_15 <- predict_glm(unisum_15$uni)
  store(5, uni_p_15, unisum_15$varsel_uni, 0.15, "UVS_15")

  xm <- as.matrix(x)
  xvm <- as.matrix(xval)
  lasso <- glmnet::cv.glmnet(xm, y, alpha = 1, family = "binomial",
                            type.measure = "deviance")
  lambda_min <- lasso$lambda.min
  lambda_1se <- lasso$lambda.1se
  lassomin_p <- as.vector(predict(lasso, xvm, s = "lambda.min", type = "response"))
  lasso1se_p <- as.vector(predict(lasso, xvm, s = "lambda.1se", type = "response"))
  varsel_min <- as.integer(as.numeric(coef(lasso, s = "lambda.min"))[-1] != 0)
  varsel_1se <- as.integer(as.numeric(coef(lasso, s = "lambda.1se"))[-1] != 0)
  store(6, lassomin_p, varsel_min, lambda_min, "LS_min")
  store(7, lasso1se_p, varsel_1se, lambda_1se, "LS_1se")

  nfolds <- 10
  f <- nfolds / (nfolds - 1) - 1
  mod_lasso <- mod_penal_ave_foreach(x = xm, y = y, bn = 20, method = "lasso",
                                    f = f, parallel = FALSE, nfolds = nfolds,
                                    boot = TRUE)
  eta_mod_lasso <- as.vector(cbind(1, xm) %*% mod_lasso$beta.boot)
  p_mod_lasso <- plogis(as.vector(cbind(1, xvm) %*% mod_lasso$beta.boot))
  varsel_mod_lasso <- as.integer(as.numeric(mod_lasso$beta.boot)[-1] != 0)
  store(8, p_mod_lasso, varsel_mod_lasso,
        as.numeric(mod_lasso$lambda.boot), "LS_mod")

  ridge <- glmnet::cv.glmnet(xm, y, alpha = 0, family = "binomial",
                            type.measure = "deviance")
  ridge_p <- as.vector(predict(ridge, xvm, s = "lambda.min", type = "response"))
  store(9, ridge_p, rep(1, n.para), ridge$lambda.min, "Ridge")

  ridge_beta <- as.numeric(coef(ridge, s = "lambda.min"))[-1]
  penalty_weights <- 1 / pmax(abs(ridge_beta), 1e-6)
  penalty_weights <- penalty_weights / mean(penalty_weights)
  alasso <- glmnet::cv.glmnet(xm, y, alpha = 1, penalty.factor = penalty_weights,
                             family = "binomial", type.measure = "deviance")
  alasso_p <- as.vector(predict(alasso, xvm, s = "lambda.min", type = "response"))
  varsel_alasso <- as.integer(as.numeric(coef(alasso, s = "lambda.min"))[-1] != 0)
  store(10, alasso_p, varsel_alasso, alasso$lambda.min, "ALASSO")

  sd_sources <- list(
    SD_MLE = eta_mle,
    SD_LS_min = as.vector(predict(lasso, xm, s = "lambda.min", type = "link")),
    SD_LS_mod = eta_mod_lasso
  )
  for (j in seq_along(sd_sources)) {
    nm <- names(sd_sources)[j]
    sd_fit <- prune_predictor_harrell(x, sd_sources[[j]], r2_threshold = 0.99)
    b <- sd_fit$coefficients
    if ("Intercept" %in% names(b)) names(b)[names(b) == "Intercept"] <- "(Intercept)"
    if (!"(Intercept)" %in% names(b)) stop("Missing step-down intercept")
    needed <- setdiff(names(b), "(Intercept)")
    if (!all(needed %in% colnames(xvm))) stop("Step-down predictors do not match")
    sd_Xval <- cbind("(Intercept)" = 1, xvm[, needed, drop = FALSE])
    sd_p <- plogis(as.vector(sd_Xval[, names(b), drop = FALSE] %*% b))
    store(10 + j, sd_p, sd_fit$kept, 0.99, nm)
  }

  prevalence_val <- mean(yval)
  nb_list[["Treat_all"]] <- data.frame(
    method = "Treat_all", threshold = thresholds,
    net_benefit = prevalence_val -
      (1 - prevalence_val) * thresholds / (1 - thresholds))
  nb_list[["Treat_none"]] <- data.frame(
    method = "Treat_none", threshold = thresholds, net_benefit = 0)

  colnames(method_result) <- c("ndev", "method", "calibration slope",
                               "calibration in the large", "auc", "Brier score",
                               paste0("varsel", seq_len(n.para)), "option")
  net_benefit_long <- do.call(rbind, nb_list)
  rownames(net_benefit_long) <- NULL
  list(method_result = method_result, net_benefit = net_benefit_long)
}
