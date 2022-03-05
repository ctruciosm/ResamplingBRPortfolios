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
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Michaud Non-Parametric Bootstrap
michaud_bootstrap <- function(x, B = 500) {
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat,  control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Conditional Factor Parametric Bootstrap
factor_parametric_bootstrap <- function(x, B = 500, n_factors = 1) {
  nobs <- nrow(x)
  p <- ncol(x)
  x_c <- scale(x, center = TRUE, scale = FALSE)[,1:p]
  eigen_decomposition <- eigen(cov(x_c))
  alpha_hat <- apply(x, 2, mean)
  betas_hat <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  factors_hat <- x_c %*% t(betas_hat)
  epsilon_hat <-  x_c - factors_hat %*% betas_hat
  Sigma_e_hat <- cov(epsilon_hat)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    epsilon_boot <- rmvnorm(nobs, rep(0,p), Sigma_e_hat)
    returns_boot <- matrix(rep(alpha_hat, nobs), ncol = p, byrow = TRUE) + factors_hat %*% betas_hat + epsilon_boot
    linear_model <- lm(as.matrix(returns_boot)~ factors_hat)
    alpha_hat_boot <- coef(linear_model)[1, ]
    beta_hat_boot <- coef(linear_model)[-1, ]
    epsilon_hat_boot <- residuals(linear_model)
    mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
    Sigma_boot_hat <- t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot)
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Conditional Factor Bootstrap
factor_bootstrap <- function(x, B = 500, n_factors = 1) {
  nobs <- nrow(x)
  p <- ncol(x)
  x_c <- scale(x, center = TRUE, scale = FALSE)[,1:p]
  eigen_decomposition <- eigen(cov(x_c))
  alpha_hat <- apply(x, 2, mean)
  betas_hat <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  factors_hat <- x_c %*% t(betas_hat)
  epsilon_hat <-  x_c - factors_hat %*% betas_hat
  Sigma_e_hat <- cov(epsilon_hat)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  
  for (j in 1:B) {
    epsilon_boot <- epsilon_hat[sample(1:nobs, nobs, replace = TRUE),]
    returns_boot <- matrix(rep(alpha_hat, nobs), ncol = p, byrow = TRUE) + factors_hat %*% betas_hat + epsilon_boot
    linear_model <- lm(as.matrix(returns_boot)~ factors_hat)
    alpha_hat_boot <- coef(linear_model)[1, ]
    beta_hat_boot <- coef(linear_model)[-1, ]
    epsilon_hat_boot <- residuals(linear_model)
    mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
    Sigma_boot_hat <- t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot)
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

