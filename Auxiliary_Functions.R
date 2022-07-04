#################################
###    Auxiliary Functions    ###
#################################
ir <- function(x) {
  mean(x)/sd(x)
}

calculate_portfolio_weights <- function(x, type = "tp", nboot, factors = NULL, riskfree = 0.005, lambda = 2) {
  
  w_estim <- markowitz_optimization(type, mu = colMeans(x), Sigma = cov(x), risk.free = riskfree, lambda)
  w_bootparam <- michaud_parametric_bootstrap(x, B = nboot, type, riskfree, lambda) 
  w_boot <- michaud_bootstrap(x, B = nboot, type, riskfree, lambda) 
  k <- POETKhat(t(x))$K1HL
  w_factor_bootparam <- factor_parametric_bootstrap(x, B = nboot, n_factors  = k, type, riskfree, lambda, factors = NULL) 
  w_factor_boot <- factor_bootstrap(x, B = nboot, n_factors  = k, type, riskfree, lambda, factors = NULL) 
  w_factor_bootparam_observed <- factor_parametric_bootstrap(x, B = nboot, n_factors  = NULL, type, riskfree, lambda, factors) 
  w_factor_boot_observed <- factor_bootstrap(x, B = nboot, n_factors  = NULL, type, riskfree, lambda, factors) 
  
  return(rbind(w_estim, w_bootparam, w_boot, w_factor_bootparam, w_factor_boot, w_factor_bootparam_observed, w_factor_boot_observed))
  
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
  #IR <- AV/SD
  SR <- (mean(x) - rf)/sd(x)
  # Adjusted Sharpe Ratio
  ASR <- SR*(1 + (moments::skewness(x)/6)*SR - ((moments::kurtosis(x) - 3)/24)*SR^2)
  # Sortino Ratio
  #SO <- AV/sqrt(mean(ifelse(x < 0, 0, x^2)))
  SO <- (mean(x) - rf)/sqrt(mean(ifelse(x - rf < 0, 0, (x - rf)^2)))
  output <- c(12*AV, sqrt(12)*SD, SR, ASR, SO)
  return(output)
}

calculate_to_sspw <- function(previous_weights, desired_weights, oos_returns, p) {
  l <- length(previous_weights)/p
  to <- rep(NA, l)
  sspw <- rep(NA, l)
  for (k in 1:l) {
    num <- previous_weights[(p*(k - 1) + 1):(p*k)]*(1 + oos_returns/100)
    den <- sum(num)
    updated_weights <- num/den
    to[k] <- sum(abs(desired_weights[(p*(k - 1) + 1):(p*k)] - updated_weights))
    sspw[k] <- sum(desired_weights[(p*(k - 1) + 1):(p*k)]^2)
  }
  return(list(to, sspw))
}