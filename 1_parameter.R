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

# Sample sizes
ndev <- 1505
ndev1 <- round(ndev/2)
ndev2 <- round(ndev/4)

# Obtain the true coefficients
opt_beta <- opt_beta(n.para, prev, c, weights1, sigma1)
beta0_1 <- opt_beta$beta0
beta_1 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights2, sigma2)
beta0_2 <- opt_beta$beta0
beta_2 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights3, sigma3)
beta0_3 <- opt_beta$beta0
beta_3 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights4, sigma4)
beta0_4 <- opt_beta$beta0
beta_4 <- opt_beta$beta1

opt_beta <- opt_beta(n.para, prev, c, weights5, sigma5)
beta0_5 <- opt_beta$beta0
beta_5 <- opt_beta$beta1
