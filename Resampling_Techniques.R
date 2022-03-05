###################################
###   Resampling Techniques     ###
###################################

library(tidyr)
library(readxl)
library(stringr)
library(dplyr)


# Michaud Parametric Bootstrap
michaud_parametric_bootstrap <- function(x, B = 500) {
  mu_hat <- apply(x, 2, mean)
  Sigma_hat <- cov(x)
  nobs <- nrow(x)
  p <- ncol(x)
  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Michaud Bootstrap
michaud_bootstrap <- function(x, B = 1000) {
  mu_hat <- apply(x, 2, mean)
  Sigma_hat <- cov(x)
  nobs <- nrow(x)
  p <- ncol(x)
  eta_hat <- scale(x, center = TRUE, scale = FALSE)[,1:p] %*% solve(chol(Sigma_hat))
  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    eta_boot <- eta_hat[sample(1:nobs,nobs, replace = TRUE),]
    returns_boot <- mu_hat + eta_boot %*% chol(Sigma_hat)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat,  control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Conditional Factor Parametric Bootstrap
factor_parametric_bootstrap <- function(x, B = 1000, n_factors = 1) {
  nobs <- nrow(x)
  p <- ncol(x)
  mu <- apply(x, 2, mean)
  x_centred <- scale(x, center = TRUE, scale = FALSE)[,1:p]
  eigen_decomposition <- eigen(cov(x_centred))
  betas <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  pca_factors <- x_centred %*% t(betas)
  epsilon <-  x_centred - pca_factors %*% betas
  Sigma_e <- cov(epsilon)
  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    epsilon_boot <- rmvnorm(nobs, rep(0,p), Sigma_e)
    returns_boot <- matrix(rep(mu, nobs), ncol = p, byrow = TRUE) + pca_factors %*% betas + epsilon_boot
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    #t(betas) %*% var(pca_factors) %*% betas  + cov(epsilon_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Conditional Factor Bootstrap
factor_bootstrap <- function(x, B = 1000, p) {
  n <- nrow(x)
  pca <- princomp(x)
  pca_factors <- pca$scores[,1]
  linear_model <- lm(as.matrix(x) ~ pca_factors)
  betas <- coef(linear_model)
  epsilon <- residuals(linear_model) 
  Sigma_e <- cov(epsilon)
  eta <- epsilon %*% solve(chol(Sigma_e))

  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- cbind(rep(1,n),pca_factors) %*% betas + eta[sample(1:n, n, replace = TRUE),] %*% chol(Sigma_e)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat,  control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}


