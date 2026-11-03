backward_pvalue <- function(data, response, p_threshold = 0.05) {
  
    # Start with the full model for backward elimination
    current_formula <- as.formula(paste(response, "~ ."))
    model <- glm(current_formula, data = data, family = binomial)
    predictors <- setdiff(names(data), response)
  
  # Main loop for stepwise selection
  while (TRUE) {
    removed <- FALSE
    
    # Backward Step: Remove predictors with p-value > threshold
      p_values <- summary(model)$coefficients[-1, 4]  # Exclude intercept
      names(p_values) <- rownames(summary(model)$coefficients)[-1]  # Get predictor names
      max_pval <- max(p_values, na.rm = TRUE)
      
      # If there is a predictor to remove with p-value > threshold
      if (max_pval > p_threshold) {
        worst_predictor <- names(which.max(p_values))
        current_formula <- as.formula(paste(response, "~", paste(setdiff(predictors, worst_predictor), collapse = " + ")))
        model <- glm(current_formula, data = data, family = binomial)
        predictors <- setdiff(predictors, worst_predictor)
        removed <- TRUE
      }
    
    
    # Stop if no further changes can be made
    if (!removed) {
      break
    }
  }
  # Get the final predictors
  all_predictors <- setdiff(names(data), response)
  back_predictors <- names(coef(model))[-1] # Exclude the intercept
  
  # Create the indicator row for variable retention
  indicator <- ifelse(all_predictors %in% back_predictors, 1, 0)
  
  # Return the final model and indicator row
  return(list(
    varsel_back = indicator,
    backmodel = model
  ))
}
