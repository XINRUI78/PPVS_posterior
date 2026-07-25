# n.true = NULL:
#   independent standard-normal predictors
# n.true supplied:
#   first n.true predictors have pairwise correlation 0.5
library(mvtnorm)
generate_ss <- function(n, n.para, beta0, beta, n.true = NULL) {
  
  # Covariance matrix
  sigma <- diag(n.para)
  
  # Introduce correlation if requested
  if (!is.null(n.true)) {
    sigma[1:n.true, 1:n.true] <- 0.5
    diag(sigma) <- 1
  }
  
  x <- rmvnorm(n, mean = rep(0, n.para), sigma = sigma)
  eta <- beta0 + x%*%beta
  p <- 1/(1+exp(-eta))
  y <- rbinom(n, 1, p)
  data <- data.frame(y,x)}

orient_draw_matrix <- function(draw_matrix, n_observations) {
  
  draw_matrix <- as.matrix(draw_matrix)
  
  if (ncol(draw_matrix) == n_observations) {
    return(draw_matrix)
  }
  
  if (nrow(draw_matrix) == n_observations) {
    return(t(draw_matrix))
  }
  
  stop(
    "Neither dimension of the prediction matrix matches ",
    "the number of validation observations."
  )
}
