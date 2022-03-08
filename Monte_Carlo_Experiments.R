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
nboot <- 30

MonteCarloExperiment <- function(MC, nassets, ntotal, nboot, option) {
  w = w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = w_comb_bootparam = w_comb_boot = w_comb_bootparam2 = matrix(NA, ncol = nassets, nrow = MC)
  sp <- setting_parameters(p = nassets) 
  
  if (option == "short_sales_")  constrains_opt <- list(type = 'minvol', constraint = 'lo')
  if (option == "bounded_")  constrains_opt <- list(type = 'minvol', LB = rep(0, nassets), UB = rep(2/nassets, nassets))
  if (option == "unc_")  constrains_opt <- list(type = 'minvol')
  
  
  for (i in 1:MC) {
    print(sprintf("Monte Carlo iterarion %d", i))
    set.seed(i + 123)
    data_sim <- dgp(nobs = ntotal, p = nassets, set_par = sp)
    returns_sim <- data_sim[[1]]
    Sigma_sim <- data_sim[[2]]
    mu_sim <- data_sim[[3]]
    # Weights
    w[i,] <- optimalPortfolio(Sigma = Sigma_sim, mu = mu_sim, control = constrains_opt) # True weights
    w_estim[i,] <- optimalPortfolio(Sigma = cov(returns_sim), mu = apply(returns_sim, 2, mean),control = constrains_opt) # Markowitz weights
    w_bootparam[i,] <- michaud_parametric_bootstrap(returns_sim, B = nboot, option_list = constrains_opt) # Michaud Parametric Bootstrap weights
    w_boot[i,] <- michaud_bootstrap(returns_sim, B = nboot, option_list = constrains_opt) # Michaud Bootstrap weights
    k <- POET::POETKhat(t(returns_sim))$K1HL
    w_factor_bootparam[i,] <- factor_parametric_bootstrap(returns_sim, B = nboot, n_factors  = k, option_list = constrains_opt) # Factor conditional Parametric Bootstrap weights
    w_factor_boot[i,] <- factor_bootstrap(returns_sim, B = nboot, n_factors  = k, option_list = constrains_opt) # Factor conditional Bootstrap weights
    #w_comb_bootparam[i,] <- combining_parametric_bootstrap(returns_sim, B = nboot, option_list = constrains_opt) # Combining Parametric Resampling
    #w_comb_bootparam2[i,] <- combining_parametric_bootstrap2(returns_sim, B = nboot, option_list = constrains_opt) # Combining Parametric Resampling
    #w_comb_boot[i,] <- combining_bootstrap(returns_sim, B = nboo, option_list = constrains_opt) # Combining Resampling
  }
  
  newfolder <- paste0("MC_d", nassets, "_n", ntotal)
  if (!dir.exists(newfolder)) {
    dir.create(file.path(getwd(), newfolder), recursive = TRUE)
  }
  
  write.csv(w, paste0(newfolder, "/", option,"w_", nassets, "_", ntotal, ".csv"))
  write.csv(w_estim, paste0(newfolder, "/", option,"w_estim", nassets, "_", ntotal, ".csv"))
  write.csv(w_bootparam, paste0(newfolder, "/", option,"w_bootparam", nassets, "_", ntotal, ".csv"))
  write.csv(w_boot, paste0(newfolder, "/", option,"w_boot", nassets, "_", ntotal, ".csv"))
  write.csv(w_factor_bootparam, paste0(newfolder, "/", option,"w_factor_bootparam", nassets, "_", ntotal, ".csv"))
  write.csv(w_factor_boot, paste0(newfolder, "/", option,"w_factor_boot", nassets, "_", ntotal, ".csv"))
  #write.csv(w_comb_bootparam, paste0(newfolder, "/", option,"w_comb_bootparam", nassets, "_", ntotal, ".csv"))
  #write.csv(w_comb_bootparam2, paste0(newfolder, "/", option,"w_comb_bootparam2", nassets, "_", ntotal, ".csv"))
  #write.csv(w_comb_boot, paste0(newfolder, "/", option,"w_comb_boot", nassets, "_", ntotal, ".csv"))
}

MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "short_sales_")
MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "bounded_")
MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "unc_")


