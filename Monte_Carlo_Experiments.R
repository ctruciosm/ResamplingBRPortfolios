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
library(Matrix)
source("DGP.R")
source("Auxiliary_Functions.R")
source("Resampling_Techniques.R")

# Setup
MC <- 1000
nassets <- 7
ntotal <- 60
nboot <- 500


MonteCarloExperiment <- function(MC, nassets, ntotal, nboot, option) {
  w = w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = w_factor_bootparam_obsfactor = w_factor_boot_obsfactor = matrix(NA, ncol = nassets, nrow = MC)
  sp <- setting_parameters(p = nassets) 
  
  if (option == "short_sales_")  constrains_opt <- list(type = 'minvol', constraint = 'lo')
  if (option == "bounded_")  constrains_opt <- list(type = 'minvol', LB = rep(0, nassets), UB = rep(2/nassets, nassets))
  if (option == "unc_")  constrains_opt <- list(type = 'minvol')
  
  SDoos_performance =  IRoos_performance = matrix(NA, ncol = 6, nrow = MC)
  for (i in 1:MC) {
    print(sprintf("Monte Carlo iterarion %d", i))
    set.seed(i + 123)
    data_sim <- dgp(nobs = ntotal, p = nassets, set_par = sp)
    returns_sim <- data_sim[[1]]
    oos_returns <- data_sim[[4]]
    Sigma_sim <- data_sim[[2]]
    mu_sim <- data_sim[[3]]
    factor_sim <- matrix(data_sim[[5]], ncol = 1)
    # Weights
    w[i,] <- optimalPortfolio(Sigma = Sigma_sim, mu = mu_sim, control = constrains_opt) # True weights
    w_estim[i,] <- optimalPortfolio(Sigma = cov(returns_sim), mu = apply(returns_sim, 2, mean),control = constrains_opt) # Markowitz weights
    w_bootparam[i,] <- michaud_parametric_bootstrap(returns_sim, B = nboot, option_list = constrains_opt) # Michaud Parametric Bootstrap weights
    w_boot[i,] <- michaud_bootstrap(returns_sim, B = nboot, option_list = constrains_opt) # Michaud Bootstrap weights
    w_factor_bootparam_obsfactor[i,] <- factor_parametric_bootstrap(returns_sim, B = nboot, n_factors  = 1, option_list = constrains_opt, factors = factor_sim) # Factor conditional Parametric Bootstrap weights
    w_factor_boot_obsfactor[i,] <- factor_bootstrap(returns_sim, B = nboot, n_factors  = 1, option_list = constrains_opt, factors = factor_sim) # Factor conditional Bootstrap weights
    
     # OoS Portfolios
    SDoos_performance[i,] <- c(apply(oos_returns %*% w[i,], 2 , sd),
                                apply(oos_returns %*% w_estim[i,], 2, sd),
                                apply(oos_returns %*% w_bootparam[i,], 2, sd),
                                apply(oos_returns %*% w_boot[i,], 2, sd),
                               apply(oos_returns %*% w_factor_bootparam_obsfactor[i,], 2, sd),
                               apply(oos_returns %*% w_factor_boot_obsfactor[i,], 2, sd))
    IRoos_performance[i,] <- c(apply(oos_returns %*% w[i,], 2 , ir),
                               apply(oos_returns %*% w_estim[i,], 2, ir),
                               apply(oos_returns %*% w_bootparam[i,], 2, ir),
                               apply(oos_returns %*% w_boot[i,], 2, ir),
                               apply(oos_returns %*% w_factor_bootparam_obsfactor[i,], 2, ir),
                               apply(oos_returns %*% w_factor_boot_obsfactor[i,], 2, ir))
 }
  colnames(SDoos_performance) <- c("True", "Estim", "ParamBoot", "Boot", "FactorParamBootObsF", "FactorBootObsF")
  colnames(IRoos_performance) <- c("True", "Estim", "ParamBoot", "Boot", "FactorParamBootObsF", "FactorBootObsF")
  
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
  write.csv(w_factor_bootparam_obsfactor, paste0(newfolder, "/", option,"w_factor_bootparam_obsfactor", nassets, "_", ntotal, ".csv"))
  write.csv(w_factor_boot_obsfactor, paste0(newfolder, "/", option,"w_factor_boot_obsfactor", nassets, "_", ntotal, ".csv"))
  write.table(SDoos_performance, paste0(newfolder, "/", option,"SDoos_performance", nassets, "_", ntotal, ".csv"), col.names = TRUE, sep = ",")
  write.table(IRoos_performance, paste0(newfolder, "/", option,"IRoos_performance", nassets, "_", ntotal, ".csv"), col.names = TRUE, sep = ",")
}

MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "short_sales_")
MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "bounded_")
#MonteCarloExperiment(MC, nassets, ntotal, nboot, option = "unc_")


