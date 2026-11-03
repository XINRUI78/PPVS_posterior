


unilogit <- function(x, y, p_threshold = 0.05) {
  
  # Ensure y is a factor for logistic regression
  #y <- as.factor(y)
  
  # Initialize storage for results
  results <- list()
  varsel_uni <- numeric(ncol(x))  # Binary row for significance
  
  # Loop through each column by number
  for (i in 1:ncol(x)) {
    
    # Fit univariate logistic regression using column index
    predictor <- x[, i]
    model <- glm(y ~ predictor, family = binomial)
    
    # Extract p-value
    summary_stats <- summary(model)$coefficients
    p_value <- summary_stats[2, 4]
    
    # Determine significance
    varsel_uni[i] <- ifelse(p_value < p_threshold, 1, 0)
  }  
  
  # Create the final model using significant variables
  significant_vars <- which(varsel_uni == 1)
  if (length(significant_vars) > 0) {
    final_formula <- as.formula(paste("y ~", paste(paste0("x[,", significant_vars, "]"), collapse = " + ")))
    final_model <- glm(final_formula, family = binomial)
  } else {
    final_model <- glm(y ~ 1, family = binomial)
  }
  
  # Return results
  return(list(
    varsel_uni = varsel_uni,
    uni = final_model
  ))
}


# Inputs

#' @param x        a matrix, the design matrix
#' @param y        a vector, the binary response variable
#' @param bn       a scalar, the number of bootstrap iterations
#' @param method  a string, values "ridge" or "lasso"Numeric where to trim (proportion)
#' @param nfolds   a scalar, the number of cross-validation folds
#' @param parallel a logical, TRUE for parallel computing
#' @param boot     a logical, TRUE means sample with replacement (the default)
#'
#' Outputs: a list with two elements.
#'          - beta.boot is the calculated vector of regression coefficients
#'          - lambda.boot is vector of beta's


mod_penal_ave_foreach <- function(x, y, bn=20, method="ridge", nfolds=10, f = nfolds/(nfolds-1)-1,
                                  parallel=TRUE,  boot=TRUE){
  library(glmnet)
  library(foreach)
  
  if (method=="lasso") a=1
  if (method=="ridge") a=0
  
  pred        <- NULL ; beta.boot   <- NULL
  lambda.boot <- NULL ; cv.boot     <- NULL
  yval        <- y
  xval        <- x
  npred       <- ncol(x)
  data        <- cbind(x,y)
  
  fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE, parallel = T)
  
  # Obtain the default sequence of lambda, and include some smaller values too
  lambdaseq <- fit$lambda
  lambdaseq <- unique(c(lambdaseq, seq( lambdaseq[length(lambdaseq)], lambdaseq[length(lambdaseq)]/20,
                                        -lambdaseq[length(lambdaseq)]/20 )))
  
  # cvd is the cross-validated deviance obtained in each of the 'over-sampled datasets'
  cvd <- NULL
  
  cv.boot <- foreach(i=1:bn, .combine='rbind', .packages=c('MASS','RcppNumerical', 'pROC', 'brglm2', 'detectseparation',
                                                           'RcppNumerical', 'logistf', 'glmnet')) %do% {
                                                             
                                                             
                                                             withwarnings <- function(expr) {
                                                               
                                                               val        <- NULL
                                                               myWarnings <- NULL
                                                               wHandler   <- function(w) {
                                                                 myWarnings <<- c(myWarnings, w$message)
                                                                 invokeRestart("muffleWarning")
                                                               }
                                                               
                                                               myError  <- NULL
                                                               eHandler <- function(e) {
                                                                 myError <<- e$message
                                                                 NULL
                                                               }
                                                               val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
                                                               list(value = val, warnings = myWarnings, error=myError)
                                                             }
                                                             
                                                             bs     <- sample(nrow(data), replace=T)
                                                             databs <- data[bs,]
                                                             
                                                             if (boot==FALSE) databs<-data
                                                             
                                                             if (f==0) databs  <- data else
                                                             {
                                                               
                                                               bs2     <- sample(nrow(data)*f, replace=T)
                                                               databs2 <- data[bs2,]
                                                               databs  <- rbind(databs,databs2)
                                                             }
                                                             
                                                             xbs   <- databs[,1:npred]
                                                             ybs   <- databs[,npred+1]
                                                             
                                                             fitbs <- withwarnings(cv.glmnet(xbs, ybs, type.measure = "deviance", family = "binomial", alpha = a, nfolds = nfolds,
                                                                                             standardize = TRUE, parallel = parallel, lambda = lambdaseq))
                                                             
                                                             
                                                             if ( length(fitbs$error) == 0 ) {
                                                               fitbs  <- fitbs$value
                                                               cvd    <- fitbs$cvm
                                                               
                                                             } else
                                                               
                                                             {
                                                               cvd <- cvd
                                                             }
                                                             
                                                             rbind(cvd)
                                                             
                                                           }
  
  b <- cbind(lambdaseq, colMeans(cv.boot, na.rm = TRUE))
  
  # Selected the value of lambda that minimises the cross-validated deviance
  
  lambda.boot <- b[order(b[,2]),] [1,1]
  
  
  # Obtain the solution using ridge/lasso at the selected value of lambda for the *original data*
  
  fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE,
                parallel = T, lambda = lambdaseq)
  
  beta.boot   <- as.vector(coef(fit, s = lambda.boot))
  
  return(list("beta.boot" = beta.boot, "lambda.boot" = lambda.boot))
}

