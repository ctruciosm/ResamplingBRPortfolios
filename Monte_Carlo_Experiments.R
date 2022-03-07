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
library(lubridate)
library(POET)
library(rrcov)                  
library(cvCovEst)
library(nlshrink)
library(Matrix)
source("DGP.R")
source("Auxiliary_Functions.R")
source("Resampling_Techniques.R")

# Setup
MC <- 1000
nassets <- 5
ntotal <- 60
nboot <- 500
w = w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = w_comb_bootparam = w_comb_boot = w_comb_bootparam2 = matrix(NA, ncol = nassets, nrow = MC)
option = "unc_mvp_"  # "short_sales_", "bounded_"
sp <- setting_parameters(p = nassets) 
for (i in 1:MC) {
  print(sprintf("Monte Carlo iterarion %d", i))
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
  w_comb_bootparam2[i,] <- combining_parametric_bootstrap2(returns_sim, B = nboot) # Combining Parametric Resampling
  w_comb_boot[i,] <- combining_bootstrap(returns_sim, B = nboot) # Combining Resampling
}


if (dir.exists(paste0("MC_d", nassets, "_n", ntotal))) {
  unlink(crytocurrency, recursive = TRUE)
  dir.create(file.path(getwd(), paste0("MC_d", nassets, "_n", ntotal)), recursive = TRUE)
} else{
  dir.create(file.path(getwd(), paste0("MC_d", nassets, "_n", ntotal)), recursive = TRUE)
}

write.csv(w, paste0(option,"w_", nassets, "_", ntotal, ".csv"))
write.csv(w_estim, paste0(option,"w_estim", nassets, "_", ntotal, ".csv"))
write.csv(w_bootparam, paste0(option,"w_bootparam", nassets, "_", ntotal, ".csv"))
write.csv(w_boot, paste0(option,"w_boot", nassets, "_", ntotal, ".csv"))
write.csv(w_factor_bootparam, paste0(option,"w_factor_bootparam", nassets, "_", ntotal, ".csv"))
write.csv(w_factor_boot, paste0(option,"w_factor_boot", nassets, "_", ntotal, ".csv"))
write.csv(w_comb_bootparam, paste0(option,"w_comb_bootparam", nassets, "_", ntotal, ".csv"))
write.csv(w_comb_bootparam2, paste0(option,"w_comb_bootparam2", nassets, "_", ntotal, ".csv"))
write.csv(w_comb_boot, paste0(option,"w_comb_boot", nassets, "_", ntotal, ".csv"))


