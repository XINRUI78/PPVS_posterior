 # Generate data
library(mvtnorm)
library(stats)
generate_ss_cor <- function(n, n.para, beta0, beta, sigma){
  x <- rmvnorm(n, mean = rep(0, n.para), sigma = sigma)
  eta <- rep(beta0, n) + x%*%beta
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
