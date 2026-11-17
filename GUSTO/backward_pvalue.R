backward_pvalue <- function(data, response, p_threshold = 0.05) {

  predictors <- setdiff(names(data), response)
  model <- glm(reformulate(predictors, response), data = data, family = binomial())

  repeat {

    if (length(predictors) == 0) break

    sm <- summary(model)$coefficients
    p_values <- sm[-1, 4]
    coef_names <- rownames(sm)[-1]

    # Match coefficient names to the original model terms
    mm <- model.matrix(model)
    assign <- attr(mm, "assign")
    term_labels <- attr(terms(model), "term.labels")

    coef_position <- match(coef_names, colnames(mm))
    term_number <- assign[coef_position]
    predictor_names <- term_labels[term_number]

    # All factors in your data are binary, so each predictor has one coefficient
    names(p_values) <- predictor_names

    p_values[is.na(p_values)] <- 1
    max_pval <- max(p_values)

    if (max_pval <= p_threshold) break

    worst_predictor <- names(which.max(p_values))

    predictors <- setdiff(predictors, worst_predictor)

    if (length(predictors) == 0) {
      model <- glm(reformulate(NULL, response), data = data, family = binomial())
    } else {
      model <- glm(reformulate(predictors, response), data = data, family = binomial())
    }
  }

  all_predictors <- setdiff(names(data), response)
  back_predictors <- predictors
  indicator <- as.integer(all_predictors %in% back_predictors)

  return(list(varsel_back = indicator, backmodel = model))
}
