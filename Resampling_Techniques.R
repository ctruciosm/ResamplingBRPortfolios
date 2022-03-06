###################################
###   Resampling Techniques     ###
###################################

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
    Sigma_boot_hat <- tryCatch(cov(returns_boot), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
      mu_boot_hat <- apply(returns_boot, 2, mean)
      Sigma_boot_hat <- tryCatch(cov(returns_boot), error = 1)
    }
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
    Sigma_boot_hat <- tryCatch(cov(returns_boot), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
      mu_boot_hat <- apply(returns_boot, 2, mean)
      Sigma_boot_hat <- tryCatch(cov(returns_boot), error = 1)
    }
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
    beta_hat_boot <- matrix(coef(linear_model)[-1, ], ncol = p)
    epsilon_hat_boot <- residuals(linear_model)
    mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
    Sigma_boot_hat <- tryCatch(t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      epsilon_boot <- rmvnorm(nobs, rep(0,p), Sigma_e_hat)
      returns_boot <- matrix(rep(alpha_hat, nobs), ncol = p, byrow = TRUE) + factors_hat %*% betas_hat + epsilon_boot
      linear_model <- lm(as.matrix(returns_boot)~ factors_hat)
      alpha_hat_boot <- coef(linear_model)[1, ]
      beta_hat_boot <- matrix(coef(linear_model)[-1, ], ncol = p)
      epsilon_hat_boot <- residuals(linear_model)
      mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
      Sigma_boot_hat <- tryCatch(t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot), error = 1)
    }
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Conditional Factor Non-Parametric Bootstrap
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
    beta_hat_boot <- matrix(coef(linear_model)[-1, ], ncol = p)
    epsilon_hat_boot <- residuals(linear_model)
    mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
    Sigma_boot_hat <- tryCatch(t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      epsilon_boot <- epsilon_hat[sample(1:nobs, nobs, replace = TRUE),]
      returns_boot <- matrix(rep(alpha_hat, nobs), ncol = p, byrow = TRUE) + factors_hat %*% betas_hat + epsilon_boot
      linear_model <- lm(as.matrix(returns_boot)~ factors_hat)
      alpha_hat_boot <- coef(linear_model)[1, ]
      beta_hat_boot <- matrix(coef(linear_model)[-1, ], ncol = p)
      epsilon_hat_boot <- residuals(linear_model)
      mu_boot_hat <- alpha_hat_boot + apply(factors_hat,2,mean) %*% beta_hat_boot 
      Sigma_boot_hat <- tryCatch(t(beta_hat_boot) %*% cov(factors_hat) %*% beta_hat_boot + cov(epsilon_hat_boot), error = 1)
    }
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Combining Parametric Bootstrap
combining_parametric_bootstrap <- function(x, B = 500) {
  mu_hat <- apply(x, 2, mean)
  Sigmas <- list(cov(x), 
                 CovMcd(x)$cov, 
                 CovMest(x)$cov,
                 CovMrcd(x)$cov,
                 CovMve(x)$cov,
                 CovOgk(x)$cov,
                 CovSest(x)$cov,
                 nlShrinkLWEst(x),  #Analytical Non-Linear Shrinkage 
                 nlshrink_cov(x),  # Non-linear Shrinkage  QUEST
                 covEstimation(x, control = list(type = 'lw')),  # Ledoit and Wolf (2003
                 covEstimation(x, control = list(type = 'oneparm')), #Ledoit and Wolf (2004)
                 covEstimation(x, control = list(type = 'const')),  # Ledoit and Wolf (2002)
                 covEstimation(x, control = list(type = 'cor')),  # Ledoit and Wolf (2003)
                 covEstimation(x, control = list(type = 'diag')),  # Ledoit and Wolf (2002)
                 covEstimation(x, control = list(type = 'large'))) # Ledoit and Wolf (2004)
  n_sigmas <- length(Sigmas)
  n_sigmas_index <- 1:n_sigmas
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    selected_sigma <- sample(n_sigmas_index,1)
    returns_boot <- rmvnorm(nobs, mu_hat, Sigmas[[selected_sigma]])
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- tryCatch(covariance_method(returns_boot, method = selected_sigma), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      returns_boot <- rmvnorm(nobs, mu_hat, Sigmas[[selected_sigma]])
      mu_boot_hat <- apply(returns_boot, 2, mean)
      Sigma_boot_hat <- tryCatch(covariance_method(returns_boot, method = selected_sigma), error = 1)
    }
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

# Combining Non-Parametric Bootstrap
combining_bootstrap <- function(x, B = 500) {
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  n_sigmas_index = 1:15  # 15 Sigma estimator implemented in covariance_method
  for (j in 1:B) {
    selected_sigma <- sample(n_sigmas_index,1)
    returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- tryCatch(covariance_method(returns_boot, method = selected_sigma), error = 1)
    while (!is.matrix(Sigma_boot_hat) | !isSymmetric.matrix(Sigma_boot_hat) | !all(eigen(Sigma_boot_hat)$values > 0)) {
      returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
      mu_boot_hat <- apply(returns_boot, 2, mean)
      Sigma_boot_hat <- tryCatch(covariance_method(returns_boot, method = selected_sigma), error = 1)
    }
    w_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(w_boot, 2, mean))
}

