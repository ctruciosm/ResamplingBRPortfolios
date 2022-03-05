###################################
###   Monte Carlo Experiments   ###
###################################

rm(list = ls())
library(mvtnorm)
library(RiskPortfolios)
library(tidyr)
library(readxl)
library(stringr)
library(dplyr)
library(POET)
library(rrcov)                  
library(cvCovEst)
library(nlshrink)
library(covmat)
library(matrixcalc)
source("DGP.R")
source("Auxiliary_Functions.R")
source("Resampling_Techniques.R")

# Setup
MC <- 10
nassets <- 10
ntotal <- 60
nboot <- 500
w = w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = w_comb_bootparam = w_comb_boot = matrix(NA, ncol = nassets, nrow = MC)

sp <- setting_parameters(p = nassets) 
for (i in 1:MC) {
  print(i)
  set.seed(i + 123)
  data_sim <- dgp(nobs = ntotal, p = nassets, set_par = sp)
  returns_sim <- data_sim[[1]]
  Sigma_sim <- data_sim[[2]]
  mu_sim <- data_sim[[3]]
  # Weights
  w[i,] <- optimalPortfolio(Sigma = Sigma_sim, mu = mu_sim, control = list(type = 'minvol')) # True weights
  w_estim[i,] <- optimalPortfolio(Sigma = cov(returns_sim), mu = apply(returns_sim, 2, mean),control = list(type = 'minvol')) # Markowitz weights
  w_bootparam[i,] <- michaud_parametric_bootstrap(returns_sim, B = nboot) # Michaud Parametric Bootstrap weights
  w_boot[i,] <- michaud_bootstrap(returns_sim, B = nboot) # Michaud Bootstrap weights
  k <- POET::POETKhat(t(returns_sim))$K1HL
  w_factor_bootparam[i,] <- factor_parametric_bootstrap(returns_sim, B = nboot, n_factors  = k) # Factor conditional Parametric Bootstrap weights
  w_factor_boot[i,] <- factor_bootstrap(returns_sim, B = nboot, n_factors  = k) # Factor conditional Bootstrap weights
  w_comb_bootparam[i,] <- combining_parametric_bootstrap(returns_sim, B = nboot) # Combining Parametric Resampling
  w_comb_boot[i,] <- combining_bootstrap(returns_sim, B = nboot) # Combining Resampling
}

write.csv(w,"w.csv")
write.csv(w_estim,"w_estim.csv")
write.csv(w_bootparam,"w_bootparam.csv")
write.csv(w_boot,"w_boot.csv")
write.csv(w_factor_bootparam,"w_factor_bootparam.csv")
write.csv(w_factor_boot,"w_factor_boot.csv")
write.csv(w_comb_bootparam,"w_comb_bootparam.csv")
write.csv(w_comb_boot,"w_comb_boot.csv")


