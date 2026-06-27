######### Set simulation parameters
n.para <- 30
prev <- 0.3
c <- 0.8
nval <- 10000
percentage1 <- c(0.1, 0.2, 0.2, 0.5)
percentage2 <- c(0.5, 0, 0, 0.5)
percentage3 <- c(0.2, 0.4, 0.4, 0)

make_weights <- function(n.para, percentage) {
  # Define number of predictors with relative strengths
  strong <- round(percentage[1] * n.para)
  medium <- round(percentage[2] * n.para)
  weak <- round(percentage[3] * n.para)
  noise <- n.para - strong - medium - weak
  
  weights <- c(
    rep(1, strong),
    rep(0.5, medium),
    rep(0.25, weak),
    rep(0, noise)
  )
  return(weights)
}


# Assign relative strengths

weights1 <- make_weights(n.para, percentage1)
weights2 <- make_weights(n.para, percentage2)
weights3 <- make_weights(n.para, percentage3)

ndev <- 1538
ndev1 <- round(ndev/2)
ndev2 <- round(ndev/4)
######### obtain the coefficents of beta
library(MASS)  # For mvrnorm to simulate predictors
library(pROC)  # For AUC calculation
opt_beta <- function(n.para, prev, c, weights) {
  # Generate predictors (X) from multivariate normal distribution
  n = 500000
  x <- mvtnorm::rmvnorm(n, mean = rep(0, n.para), sigma = diag(n.para))
  
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
  initial_para <- c(-2, 1)
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


opt_beta <- opt_beta(n.para, prev, c, weights1)
beta0_1 <- opt_beta$beta0
beta_1 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights2)
beta0_2 <- opt_beta$beta0
beta_2 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights3)
beta0_3 <- opt_beta$beta0
beta_3 <- opt_beta$beta1

# Define the optimizer function
opt_beta_var <- function(n.para, n.true, prev, c, weights) {
  # Generate predictors (X) from multivariate normal distribution
  n = 500000
  sigma <- diag(n.para)
  sigma[1:n.true, 1:n.true] <- 0.5
  diag(sigma) <- 1  # Ensure variances remain 1
  
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

n.true <- 15
opt_beta <- opt_beta_var(n.para, n.true, prev, c, weights1)
beta0_4 <- opt_beta$beta0
beta_4 <- opt_beta$beta1