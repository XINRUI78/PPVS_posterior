
unilogit <- function(x, y, p_threshold = 0.05) {

  x <- as.data.frame(x)
  predictors <- names(x)
  dat <- data.frame(y = y, x)

  varsel_uni <- integer(ncol(x))

  for (i in seq_along(predictors)) {
    fm <- reformulate(predictors[i], response = "y")
    model <- glm(fm, data = dat, family = binomial())

    p_value <- summary(model)$coefficients[2, 4]
    varsel_uni[i] <- as.integer(!is.na(p_value) && p_value < p_threshold)
  }

  selected <- predictors[varsel_uni == 1]

  if (length(selected)) {
    final_formula <- reformulate(selected, response = "y")
  } else {
    final_formula <- y ~ 1
  }

  final_model <- glm(final_formula, data = dat, family = binomial())

  list(varsel_uni = varsel_uni, uni = final_model)
}
