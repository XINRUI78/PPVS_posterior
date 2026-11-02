library(mvtnorm)
library(pROC)
opt_beta <- function(n.para, prev, c, weights, sigma) {
  # Generate predictors (X) from multivariate normal distribution
  n = 500000
  # Generate data
  x <- mvtnorm::rmvnorm(n, mean = rep(0, n.para), sigma = sigma)
  
  
  objective <- function(para){
    beta0 <- para[1]  # Intercept
    s <- para[2]      # Scaling factor
    beta1 <- s * weights
    eta <- rep(beta0, n) + x %*% beta1
    p <- 1/(1+exp(-eta))
    y <- stats::rbinom(n, 1, p)
    pest <- mean(y)
    cstat <- pROC::roc(response = as.vector(y), predictor = as.vector(p), levels = c(0, 1), direction = "<")
    cest <- as.vector(cstat$auc)
    return((pest - prev)^2 + (cest - c)^2)
  }
  # Initial guesses for beta0 and s
  initial_para <- c(log(prev/(1-prev)), 1)
  tol = 1e-6
  # Perform optimization
  result <- optim(
    par = initial_para,
    fn = objective,
    method = "Nelder-Mead", 
    control = list(abstol = tol)
  )
  
  # Extract optimized coefficients
  beta0_opt <- result$par[1]
  s_opt <- result$par[2]
  beta1_opt <- s_opt * weights
  
  list(
    beta0 = beta0_opt,
    beta1 = beta1_opt,
    s = s_opt
  )
}
