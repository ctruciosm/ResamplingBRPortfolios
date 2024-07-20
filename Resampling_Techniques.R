###################################
###   Resampling Techniques     ###
###################################

# Michaud Parametric Bootstrap
michaud_parametric_bootstrap <- function(x, B = 500, type = "tp", riskfree = 0.005, lambda = 2) {
  mu_hat <- colMeans(x)
  Sigma_hat <- cov(x)
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
    mu_boot_hat <- colMeans(returns_boot)
    c_returns_boot <- scale(returns_boot, center = T, scale = F)
    Sigma_boot_hat <- cov(c_returns_boot)
    Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
    w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) || any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      returns_boot <- rmvnorm(nobs, mu_hat, Sigma_hat)
      mu_boot_hat <- colMeans(returns_boot)
      c_returns_boot <- scale(returns_boot, center = T, scale = F)
      Sigma_boot_hat <- cov(c_returns_boot)
      Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
      w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  
  return(list(w = round(colMeans(w_boot),6), w_sd = apply(w_boot, 2, sd)))
}

# Michaud Non-Parametric Bootstrap
michaud_bootstrap <- function(x, B = 500, type = "tp", riskfree = 0.005, lambda = 2) {
  nobs <- nrow(x)
  p <- ncol(x)
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    returns_boot <- x[sample(1:nobs, nobs, replace = TRUE),]
    mu_boot_hat <- colMeans(returns_boot)
    c_returns_boot <- scale(returns_boot, center = T, scale = F)
    Sigma_boot_hat <- cov(c_returns_boot)
    Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
    w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) || any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      returns_boot <- x[sample(1:nobs,nobs, replace = TRUE),]
      mu_boot_hat <- colMeans(returns_boot)
      c_returns_boot <- scale(returns_boot, center = T, scale = F)
      Sigma_boot_hat <- cov(c_returns_boot)
      Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
      w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }

  return(list(w = round(colMeans(w_boot),6), w_sd = apply(w_boot, 2, sd)))
}

# Conditional Factor Parametric Bootstrap
factor_parametric_bootstrap <- function(x, B = 500, n_factors = 1, type = "tp", riskfree = 0.005, lambda = 2, factors = NULL) {
  nobs <- nrow(x)
  p <- ncol(x)
  if (is.null(factors)) {
    x_c <- scale(x, center = TRUE, scale = FALSE)
    eigen_decomposition <- eigen(cov(x_c))
    eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
    factors <- x %*% t(eigen_vectors)
  }
  n_factors <- ncol(factors)
  
  model <- lm(as.matrix(x)[-1, ] ~ as.matrix(factors[1:(nobs - 1), 1:n_factors]))
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_e_hat <- (1/model$df.residual) * (nobs - 2) * cov(epsilon_hat)
  Sigma_F <- ifelse(n_factors == 1, var(factors[1:(nobs - 1), ]), cov(matrix(factors[1:(nobs - 1), ], ncol = n_factors)))
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    epsilon_boot <- rmvnorm(nobs - 1, rep(0,p), Sigma_e_hat)
    returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors[1:(nobs - 1), ] %*% beta_hat + epsilon_boot
    model_boot <- lm(as.matrix(returns_boot)~ factors[1:(nobs - 1), ])
    alpha_hat_boot <- coef(model_boot)[1, ]
    beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
    epsilon_hat_boot <-  model_boot$residuals
    mu_boot_hat <- as.numeric(alpha_hat_boot + apply(matrix(factors[1:(nobs - 1), ], ncol = n_factors), 2, mean) %*% beta_hat_boot) 
    Sigma_e_boot_hat <- (1/model_boot$df.residual) * (nobs - 2) * cov(epsilon_hat_boot)
    Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
    Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
    w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    
    while (any(is.na(w_boot[j, ])) || any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      epsilon_boot <- rmvnorm(nobs - 1, rep(0,p), Sigma_e_hat)
      returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors[1:(nobs - 1), ] %*% beta_hat + epsilon_boot
      model_boot <- lm(as.matrix(returns_boot)~ factors[1:(nobs - 1), ])
      alpha_hat_boot <- coef(model_boot)[1, ]
      beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
      epsilon_hat_boot <-  model_boot$residuals
      mu_boot_hat <- as.numeric(alpha_hat_boot + apply(matrix(factors[1:(nobs - 1), ], ncol = n_factors), 2, mean) %*% beta_hat_boot) 
      Sigma_e_boot_hat <- (1/model_boot$df.residual) * (nobs - 2) * cov(epsilon_hat_boot)
      Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
      Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
      w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  
  return(list(w = round(colMeans(w_boot),6), w_sd = apply(w_boot, 2, sd)))
}

# Conditional Factor Non-Parametric Bootstrap
factor_bootstrap <- function(x, B = 500, n_factors = 1, type = "tp", riskfree = 0.005, lambda = 2, factors = NULL) {
  nobs <- nrow(x)
  p <- ncol(x)
  if (is.null(factors)) {
    x_c <- scale(x, center = TRUE, scale = FALSE)
    eigen_decomposition <- eigen(cov(x_c))
    eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
    factors <- x %*% t(eigen_vectors)
  } else {
    n_factors <- ncol(factors)
  }
  model <- lm(as.matrix(x)[-1, ] ~ factors[1:(nobs - 1), ])
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_F <- ifelse(n_factors == 1, var(factors[1:(nobs - 1), ]), cov(matrix(factors[1:(nobs - 1), ], ncol = n_factors)))
  w_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    epsilon_boot <- epsilon_hat[sample(1:(nobs - 1), nobs - 1, replace = TRUE),]
    returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors[1:(nobs - 1), ] %*% beta_hat + epsilon_boot
    model_boot <- lm(as.matrix(returns_boot)~ factors[1:(nobs - 1), ])
    alpha_hat_boot <- coef(model_boot)[1, ]
    beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
    epsilon_hat_boot <- model_boot$residuals
    mu_boot_hat <- as.numeric(alpha_hat_boot + apply(matrix(factors[1:(nobs - 1), ], ncol = n_factors), 2, mean) %*% beta_hat_boot) 
    Sigma_e_boot_hat <- (1/model_boot$df.residual) * (nobs - 2) * cov(epsilon_hat_boot)
    Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
    Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
    w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    while (any(is.na(w_boot[j, ])) || any(eigen(Sigma_boot_hat)$values <= 1e-08)) {
      epsilon_boot <- epsilon_hat[sample(1:(nobs - 1), nobs - 1, replace = TRUE),]
      returns_boot <- matrix(rep(alpha_hat, nobs - 1), ncol = p, byrow = TRUE) + factors[1:(nobs - 1), ] %*% beta_hat + epsilon_boot
      model_boot <- lm(as.matrix(returns_boot)~ factors[1:(nobs - 1), ])
      alpha_hat_boot <- coef(model_boot)[1, ]
      beta_hat_boot <- matrix(coef(model_boot)[-1, ], ncol = p)
      epsilon_hat_boot <- model_boot$residuals
      mu_boot_hat <- as.numeric(alpha_hat_boot + apply(matrix(factors[1:(nobs - 1), ], ncol = n_factors), 2, mean) %*% beta_hat_boot) 
      Sigma_e_boot_hat <- (1/model_boot$df.residual) * (nobs - 2) * cov(epsilon_hat_boot)
      Sigma_boot_hat <- t(beta_hat_boot) %*% Sigma_F %*% beta_hat_boot + Sigma_e_boot_hat
      Sigma_boot_hat <- 0.5*(Sigma_boot_hat + t(Sigma_boot_hat))
      w_boot[j, ] <- tryCatch({markowitz_optimization(type, mu = mu_boot_hat, Sigma = Sigma_boot_hat, risk.free = riskfree, lambda)}, warning = function(w) rep(NA, p), error = function(e) rep(NA, p)) 
    }
  }
  
  return(list(w = round(colMeans(w_boot),6), w_sd = apply(w_boot, 2, sd)))
}

