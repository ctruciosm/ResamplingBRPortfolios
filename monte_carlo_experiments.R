###################################
###   Monte Carlo Experiments   ###
###################################

rm(list = ls())
library(mvtnorm)
library(RiskPortfolios)
library(POET)  # Bai and NG  / Hallin and Liska
source("auxiliary_functions.R")

# Setup
MC <- 1000
nassets <- 10
ntotal <- 6000
nboot <- 1000
w = w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = matrix(NA, ncol = nassets, nrow = MC)

sp <- setting_parameters(p = nassets) 
for (i in 1:MC) {
  set.seed(i + 123)
  data_sim <- dgp(n = ntotal, p = nassets, set_par = sp)
  returns_sim <- data_sim[[1]]
  Sigma_sim <- data_sim[[2]]
  mu_sim <- data_sim[[3]]
  # Weights
  w[i,] <- optimalPortfolio(Sigma = Sigma_sim, mu = mu_sim, control = list(type = 'minvol')) # True weights
  w_estim[i,] <- optimalPortfolio(Sigma = cov(returns_sim), mu = apply(returns_sim, 2, mean),control = list(type = 'minvol')) # Markowitz weights
  w_bootparam[i,] <- michaud_parametric_bootstrap(returns_sim, B = nboot, p = nassets) # Michaud Parametric Bootstrap weights
  w_boot[i,] <- michaud_bootstrap(returns_sim, B = nboot, p = nassets) # Michaud Bootstrap weights
  w_factor_bootparam[i,] <- factor_parametric_bootstrap(returns_sim, B = nboot, p = nassets) # Factor conditional Parametric Bootstrap weights
  w_factor_boot[i,] <- factor_bootstrap(returns_sim, B = nboot, p = nassets) # Factor conditional Bootstrap weights
}

write.csv(w,"w.csv")
write.csv(w_estim,"w_estim.csv")
write.csv(w_bootparam,"w_bootparam.csv")
write.csv(w_boot,"w_boot.csv")
write.csv(w_factor_bootparam,"w_factor_bootparam.csv")
write.csv(w_factor_boot,"w_factor_boot.csv")


