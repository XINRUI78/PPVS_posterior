source("opt_beta.R")
######### Set simulation parameters
n.para <- 30
prev <- 0.3
c <- 0.8
nval <- 10000

percentage1 <- c(0.1, 0.2, 0.2, 0.5)
percentage2 <- c(0.1, 0.2, 0.2, 0.5)
percentage3 <- c(0.5, 0, 0, 0.5)
percentage4 <- c(0.2, 0.4, 0.4, 0)
percentage5 <- c(0.1, 0.2, 0.2, 0.5)

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
weights4 <- make_weights(n.para, percentage4)
weights5 <- make_weights(n.para, percentage5)

# Assign correlation matrix
sigma1 <- diag(n.para)
sigma2 <- diag(n.para)
sigma2[1:n.para, 1:n.para] <- 0.5
diag(sigma2) <- 1  # Ensure variances remain 1
sigma3 <- diag(n.para)
sigma4 <- diag(n.para)
sigma5 <- diag(n.para)
for (j in seq(1, 7, by=2)) {
       sigma5[j, j+1] <- -0.7
       sigma5[j+1, j] <- -0.7
}

sigma_list <- c(sigma1, sigma2, sigma3, sigma4, sigma5)

# Sample sizes
ndev <- 1505
ndev1 <- round(ndev/2)
ndev2 <- round(ndev/4)

# Obtain the true coefficients

# Obtain the true coefficients

result1 <- opt_beta(n.para, prev, c, weights1, sigma1)
beta0_1 <- result1$beta0
beta_1  <- result1$beta1

result2 <- opt_beta(n.para, prev, c, weights2, sigma2)
beta0_2 <- result2$beta0
beta_2  <- result2$beta1

result3 <- opt_beta(n.para, prev, c, weights3, sigma3)
beta0_3 <- result3$beta0
beta_3  <- result3$beta1

result4 <- opt_beta(n.para, prev, c, weights4, sigma4)
beta0_4 <- result4$beta0
beta_4  <- result4$beta1

result5 <- opt_beta(n.para, prev, c, weights5, sigma5)
beta0_5 <- result5$beta0
beta_5  <- result5$beta1

beta_list <- list(beta_1, beta_2, beta_3, beta_4, beta_5)
beta0_list <- list(beta0_1, beta0_2, beta0_3, beta0_4, beta0_5)
