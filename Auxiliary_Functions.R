#################################
###    Auxiliary Functions    ###
#################################
ir <- function(x) {
  mean(x)/sd(x)
}


calculate_portfolio_weights <- function(x, type = "tp", nboot, factors = NULL, riskfree = 0.005, lambda = 2) {
  nobs <- nrow(x)
  p <- ncol(x)
  
  w_estim <- markowitz_optimization(type, mu = colMeans(x), Sigma = cov_estim(x), risk.free = riskfree, lambda)
  w_bootparam <- michaud_parametric_bootstrap(x, B = nboot, type, riskfree, lambda)
  w_boot <- michaud_bootstrap(x, B = nboot, type, riskfree, lambda) 
  
  k <- POETKhat(t(x))$K1HL
  x_c <- scale(x, center = TRUE, scale = FALSE)
  eigen_decomposition <- eigen(cov_estim(x_c))
  eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:k]), nrow = k)
  pca_factors <- x %*% t(eigen_vectors)

  model <- lm(as.matrix(x)[-1, ] ~ head(pca_factors, nobs - 1))
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_e_hat <- (1/model$df.residual) * (nobs - 2) * cov_estim(epsilon_hat)
  Sigma_F <- ifelse(k == 1, var(pca_factors[1:(nobs - 1)]), cov_estim(matrix(head(pca_factors, nobs - 1), ncol = k)))
  Sigma_hat <- t(beta_hat) %*% Sigma_F %*% beta_hat + Sigma_e_hat
  mu_hat <- as.numeric(alpha_hat + apply(matrix(head(pca_factors, nobs - 1), ncol = k), 2, mean) %*% beta_hat) 
  w_estim_pca <- markowitz_optimization(type, mu = mu_hat, Sigma = Sigma_hat, risk.free = riskfree, lambda)
  w_factor_bootparam <- factor_parametric_bootstrap(x, B = nboot, n_factors  = k, type, riskfree, lambda, factors = NULL) 
  w_factor_boot <- factor_bootstrap(x, B = nboot, n_factors  = k, type, riskfree, lambda, factors = NULL) 
  
  k <- ncol(factors)
  model <- lm(as.matrix(x)[-1, ] ~ head(factors, nobs - 1))
  alpha_hat <- coef(model)[1,]
  beta_hat <- matrix(coef(model)[-1, ], ncol = p)
  epsilon_hat <-  model$residuals
  Sigma_e_hat <- (1/model$df.residual) * (nobs - 2) * cov_estim(epsilon_hat)
  Sigma_F <- ifelse(k == 1, var(factors[1:(nobs - 1)]), cov_estim(matrix(head(factors, nobs - 1), ncol = k)))
  Sigma_hat <- t(beta_hat) %*% Sigma_F %*% beta_hat + Sigma_e_hat
  mu_hat <- as.numeric(alpha_hat + apply(matrix(head(factors, nobs - 1), ncol = k), 2, mean) %*% beta_hat) 
  w_estim_observed <- markowitz_optimization(type, mu = mu_hat, Sigma = Sigma_hat, risk.free = riskfree, lambda)
  w_factor_bootparam_observed <- factor_parametric_bootstrap(x, B = nboot, n_factors  = NULL, type, riskfree, lambda, factors) 
  w_factor_boot_observed <- factor_bootstrap(x, B = nboot, n_factors  = NULL, type, riskfree, lambda, factors) 
  
  return(weights = rbind(w_estim, w_bootparam$w, w_boot$w, w_estim_pca, w_factor_bootparam$w, w_factor_boot$w, w_estim_observed, w_factor_bootparam_observed$w, w_factor_boot_observed$w))
  
}

markowitz_optimization <- function(type = "mv", mu = NULL, Sigma, risk.free = 0.05, lambda = 2) {
  if (type == "mv") weights = minvar_portfolio(Sigma)
  if (type == "tp") weights = tangency_portfolio(mu, Sigma, risk.free)
  if (type == "ef") weights = ef_portfolio(mu, Sigma, lambda)
  return(weights)
}

tangency_portfolio <- function(mu, cov.mat, risk.free) {
  # if short-selling are allowed:
  #  weights <- solve(cov.mat) %*% (er - risk.free) 
  #  weights <- as.vector(weights/sum(weights))
  #Dmat <- cov.mat
  #N <- nrow(Dmat)
  #dvec <- rep(0, N)
  #excess <- mu - risk.free
  #Amat <- cbind(-excess, diag(1,N))
  #meq <- 1 
  #bvec <- rbind(-1, matrix(0, N, 1))
  #weights <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  #weights <- weights/sum(weights)

  Dmat <- 2*cov.mat
  N <- nrow(Dmat)
  dvec <- rep(0, N)
  excess <- mu - risk.free
  Amat <- cbind(excess, diag(1,N))
  meq <- 1 
  bvec <- rbind(1, matrix(0, N, 1))
  weights <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  weights <- weights/sum(weights)
  
  return(round(weights,6))
}

ef_portfolio <- function(mu, cov.mat, lambda) {
  N <- ncol(cov.mat)
  Dmat <- lambda*cov.mat
  dvec <- mu
  Amat <- cbind(rep(1,N), diag(1,N))
  meq <- 1 
  bvec <- rbind(1, matrix(0, N, 1)) 
  weights <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  return(round(weights,6))
}

minvar_portfolio <- function(cov.mat) {
  N <- ncol(cov.mat)
  Dmat <- cov.mat
  dvec <- rep(0, N)
  Amat <- cbind(rep(1,N), diag(1,N))
  meq <- 1 
  bvec <- rbind(1, matrix(0, N, 1)) 
  weights <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  return(round(weights,6))
}

medidas <- function(x, rf = 0.5) {
  # Annualized Average
  AV <- mean(x)
  # Annualized SD
  SD <- sd(x)
  # Information (or Sharpe) Ratio
  SR <- (mean(x) - rf)/sd(x)
  # Adjusted Sharpe Ratio
  ASR <- SR*(1 + (moments::skewness(x)/6)*SR - ((moments::kurtosis(x) - 3)/24)*SR^2)
  # Sortino Ratio
  SO <- (mean(x) - rf)/sqrt(mean(ifelse(x - rf < 0, 0, (x - rf)^2)))
  output <- c(12*AV, sqrt(12)*SD, sqrt(12)*SR, sqrt(12)*ASR, sqrt(12)*SO)
  return(output)
}

calculate_to<- function(previous_weights, desired_weights, oos_returns, p) {
  oos_returns[is.na(oos_returns)] <- 0
  num <- previous_weights*(1 + oos_returns/100)
  den <- sum(num, na.rm = TRUE)
  updated_weights <- num/den
  to <- sum(abs(desired_weights - updated_weights))
  return(to)
}