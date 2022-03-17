###################################
###   Resampling Techniques     ###
###################################

# Michaud Parametric Bootstrap
michaud_parametric_bootstrap <- function(x, B = 500, option_list = list(type = 'minvol')) {
  mu_hat <- colMeans(x)
  Sigma_hat <- cov(x)
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
    mu_boot_hat <- colMeans(returns_boot)
    c_returns_boot <- scale(returns_boot, center = T, scale = F)
    Sigma_boot_hat <- 1/(nobs - 1) * t(c_returns_boot) %*% c_returns_boot
    w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) | !isSymmetric.matrix(Sigma_boot_hat) | any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
      mu_boot_hat <- colMeans(returns_boot)
      c_returns_boot <- scale(returns_boot, center = T, scale = F)
      Sigma_boot_hat <- 1/(nobs - 1) * t(c_returns_boot) %*% c_returns_boot
      w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  return(w = colMeans(w_boot))
}

# Michaud Non-Parametric Bootstrap
michaud_bootstrap <- function(x, B = 500, option_list = list(type = 'minvol')) {
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- x[sample(1:nobs, nobs, replace = TRUE),]
    mu_boot_hat <- colMeans(returns_boot)
    c_returns_boot <- scale(returns_boot, center = T, scale = F)
    Sigma_boot_hat <- 1/(nobs - 1) * t(c_returns_boot) %*% c_returns_boot
    w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) | !isSymmetric.matrix(Sigma_boot_hat) | any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
      mu_boot_hat <- colMeans(returns_boot)
      c_returns_boot <- scale(returns_boot, center = T, scale = F)
      Sigma_boot_hat <- 1/(nobs - 1) * t(c_returns_boot) %*% c_returns_boot
      w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  return(w = colMeans(w_boot))
}

# Conditional Factor Parametric Bootstrap
factor_parametric_bootstrap <- function(x, B = 500, n_factors = 1, option_list = list(type = 'minvol')) {
  nobs <- nrow(x)
  p <- ncol(x)
  x_c <- scale(x, center = TRUE, scale = FALSE)
  eigen_decomposition <- eigen(1/(nobs - 1) * t(x_c) %*% x_c)
  eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  factors_hat <- x %*% t(eigen_vectors)
  model <- lm(as.matrix(x[-1, ]) ~ factors_hat[-nobs, ])
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_e_hat <- (1/model$df.residual) * t(epsilon_hat) %*% epsilon_hat
  Sigma_F <- cov(matrix(factors_hat[-nobs, ], ncol = n_factors))
  
  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    epsilon_boot <- rmvnorm(nobs - 1, rep(0,p), Sigma_e_hat)
    returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors_hat[-nobs, ] %*% beta_hat + epsilon_boot
    model_boot <- lm(as.matrix(returns_boot)~ factors_hat[-nobs, ])
    alpha_hat_boot <- coef(model_boot)[1, ]
    beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
    epsilon_hat_boot <-  model_boot$residuals
    mu_boot_hat <- alpha_hat_boot + apply(matrix(factors_hat[-nobs, ], ncol = n_factors),2,mean) %*% beta_hat_boot 
    Sigma_e_boot_hat <- (1/model_boot$df.residual) * t(epsilon_hat_boot) %*% epsilon_hat_boot
    Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
    w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) | !isSymmetric.matrix(Sigma_boot_hat) | any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      epsilon_boot <- rmvnorm(nobs - 1, rep(0,p), Sigma_e_hat)
      returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors_hat[-nobs, ] %*% beta_hat + epsilon_boot
      model_boot <- lm(as.matrix(returns_boot)~ factors_hat[-nobs, ])
      alpha_hat_boot <- coef(model_boot)[1, ]
      beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
      epsilon_hat_boot <-  model_boot$residuals
      mu_boot_hat <- alpha_hat_boot + apply(matrix(factors_hat[-nobs, ], ncol = n_factors),2,mean) %*% beta_hat_boot 
      Sigma_e_boot_hat <- (1/model_boot$df.residual) * t(epsilon_hat_boot) %*% epsilon_hat_boot
      Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
      #Sigma_boot_hat <- tryCatch({t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat}, warning = function(w) 1, error = function(e) 1)
      w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  return(w = colMeans(w_boot))
}

# Conditional Factor Non-Parametric Bootstrap
factor_bootstrap <- function(x, B = 500, n_factors = 1, option_list = list(type = 'minvol')) {
  nobs <- nrow(x)
  p <- ncol(x)
  x_c <- scale(x, center = TRUE, scale = FALSE)
  eigen_decomposition <- eigen(1/(nobs - 1) * t(x_c) %*% x_c)
  eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  factors_hat <- x %*% t(eigen_vectors)
  model <- lm(as.matrix(x[-1, ]) ~ factors_hat[-nobs, ])
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_F <- cov(matrix(factors_hat[-nobs, ], ncol = n_factors))

  w_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    epsilon_boot <- epsilon_hat[sample(1:(nobs - 1), nobs - 1, replace = TRUE),]
    returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors_hat[-nobs, ] %*% beta_hat + epsilon_boot
    model_boot <- lm(as.matrix(returns_boot)~ factors_hat[-nobs, ])
    alpha_hat_boot <- coef(model_boot)[1, ]
    beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
    epsilon_hat_boot <- model_boot$residuals
    mu_boot_hat <- alpha_hat_boot + apply(matrix(factors_hat[-nobs, ], ncol = n_factors),2,mean) %*% beta_hat_boot 
    Sigma_e_boot_hat <- (1/model_boot$df.residual) * t(epsilon_hat_boot) %*% epsilon_hat_boot
    Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
    w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) | !isSymmetric.matrix(Sigma_boot_hat) | any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      epsilon_boot <- epsilon_hat[sample(1:(nobs - 1), nobs - 1, replace = TRUE),]
      returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors_hat[-nobs, ] %*% beta_hat + epsilon_boot
      model_boot <- lm(as.matrix(returns_boot)~ factors_hat[-nobs, ])
      alpha_hat_boot <- coef(model_boot)[1, ]
      beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
      epsilon_hat_boot <- model_boot$residuals
      mu_boot_hat <- alpha_hat_boot + apply(matrix(factors_hat[-nobs, ], ncol = n_factors),2,mean) %*% beta_hat_boot 
      Sigma_e_boot_hat <- (1/model_boot$df.residual) * t(epsilon_hat_boot) %*% epsilon_hat_boot
      Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
      w_boot[j, ] <- tryCatch({optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = option_list)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  return(w = colMeans(w_boot))
}

