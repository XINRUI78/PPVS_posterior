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
