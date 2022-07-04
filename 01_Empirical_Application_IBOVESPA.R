#################################
###   Empirical Application   ###
#################################
setwd("/Volumes/CTRUCIOS_SD/Ongoing Research/Resampling Portfolios/ResamplingBRPortfolios")
rm(list = ls())
library(mvtnorm)
library(POET)
library(Matrix)
library(kableExtra)
source("Resampling_Techniques.R")
source("Auxiliary_Functions.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#            Data              #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
monthly_data <- read.csv("monthly_data.csv")
ibovespa <- read.csv("ibovespa.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#       General Settings       #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
InS <- 120
OoS <- nrow(monthly_data) - InS
p <- ncol(monthly_data) 
nboot <- 1000
nmethods <- 7
rf <- 0.5
lambda <- 2
nconstrains <- 3

w_estim = w_bootparam = w_boot = w_factor_bootparam_pca = w_factor_boot_pca = w_factor_bootparam_obs = w_factor_boot_obs = matrix(NA, nrow = OoS, ncol = nconstrains*p)
Rport <- matrix(NA, nrow = OoS, ncol = nconstrains*nmethods)
to <- matrix(0, nrow = OoS - 1, ncol = nconstrains*nmethods)
sspw <- matrix(0, nrow = OoS, ncol = nconstrains*nmethods)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#    Out-of-sample exercise    #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
for (i in 1:OoS) { 
  print(i)
  
  r_ins <- as.matrix(monthly_data[i:(InS - 1 + i), ])
  ibov <- as.matrix(ibovespa[i:(InS - 1 + i), ])
  r_oos <- as.numeric(monthly_data[InS - 1 + i + 1, ]) 
  
  set.seed(i + 1234)
  w_mv <- calculate_portfolio_weights(x = r_ins, type = "mv",  nboot, factors = ibov)
  w_tp <- calculate_portfolio_weights(x = r_ins, type = "tp",  nboot, factors = ibov, riskfree = rf)
  w_ef <- calculate_portfolio_weights(x = r_ins, type = "ef",  nboot, factors = ibov, lambda = lambda)
  
  
  w_estim[i, ] <- c(w_mv[1, ],  w_tp[1, ], w_ef[1, ])
  w_bootparam[i, ] <- c(w_mv[2, ], w_tp[2, ], w_ef[2, ])
  w_boot[i, ] <- c(w_mv[3, ], w_tp[3, ], w_ef[3, ])
  w_factor_bootparam_pca[i, ] <- c(w_mv[4, ], w_tp[4, ], w_ef[5, ])
  w_factor_boot_pca[i, ] <- c(w_mv[5, ], w_tp[5, ], w_ef[5, ])
  w_factor_bootparam_obs[i, ] <- c(w_mv[6, ], w_tp[6, ], w_ef[6, ])
  w_factor_boot_obs[i, ] <- c(w_mv[7, ], w_tp[7, ], w_ef[7, ])
  
  if (i == 1) {
    sspw[1, ] <- c(sum(w_mv[1, ]^2),  sum(w_tp[1, ]^2), sum(w_ef[1, ]^2),
                   sum(w_mv[2, ]^2),  sum(w_tp[2, ]^2), sum(w_ef[2, ]^2),
                   sum(w_mv[3, ]^2),  sum(w_tp[3, ]^2), sum(w_ef[3, ]^2),
                   sum(w_mv[4, ]^2),  sum(w_tp[4, ]^2), sum(w_ef[4, ]^2),
                   sum(w_mv[5, ]^2),  sum(w_tp[5, ]^2), sum(w_ef[5, ]^2),
                   sum(w_mv[6, ]^2),  sum(w_tp[6, ]^2), sum(w_ef[6, ]^2),
                   sum(w_mv[7, ]^2),  sum(w_tp[7, ]^2), sum(w_ef[7, ]^2))
  } else {
    w_measures_estim <- calculate_to_sspw(w_estim[i - 1,], w_estim[i,], r_oos, p)
    w_measures_bootparam <- calculate_to_sspw(w_bootparam[i - 1,], w_bootparam[i,], r_oos, p)
    w_measures_boot <- calculate_to_sspw(w_boot[i - 1,], w_boot[i,], r_oos, p)
    w_measures_factor_bootparam_pca <- calculate_to_sspw(w_factor_bootparam_pca[i - 1,], w_factor_bootparam_pca[i,], r_oos, p)
    w_measures_factor_boot_pca <- calculate_to_sspw(w_factor_boot_pca[i - 1,], w_factor_boot_pca[i,], r_oos, p)
    w_measures_factor_bootparam_obs <- calculate_to_sspw(w_factor_bootparam_obs[i - 1,], w_factor_bootparam_obs[i,], r_oos, p)
    w_measures_factor_boot_obs <- calculate_to_sspw(w_factor_boot_obs[i - 1,], w_factor_boot_obs[i,], r_oos, p)
    to[i - 1,] <- c(w_measures_estim[[1]], w_measures_bootparam[[1]], w_measures_boot[[1]], w_measures_factor_bootparam_pca[[1]], w_measures_factor_boot_pca[[1]], w_measures_factor_bootparam_obs[[1]], w_measures_factor_boot_obs[[1]])
    sspw[i,] <- c(w_measures_estim[[2]], w_measures_bootparam[[2]], w_measures_boot[[2]], w_measures_factor_bootparam_pca[[2]], w_measures_factor_boot_pca[[2]], w_measures_factor_bootparam_obs[[2]], w_measures_factor_boot_obs[[2]])
  }
  
  Rport[i,] <-  c(r_oos %*% w_estim[i, 1:p], r_oos %*% w_bootparam[i, 1:p], r_oos %*% w_boot[i, 1:p], r_oos %*% w_factor_bootparam_pca[i, 1:p], r_oos %*% w_factor_boot_pca[i, 1:p], r_oos %*% w_factor_bootparam_obs[i, 1:p], r_oos %*% w_factor_boot_obs[i, 1:p],
                  r_oos %*% w_estim[i, (p + 1):(2 * p)], r_oos %*% w_bootparam[i, (p + 1):(2 * p)], r_oos %*% w_boot[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam_pca[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot_pca[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam_obs[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot_obs[i, (p + 1):(2 * p)],
                  r_oos %*% w_estim[i, (2*p + 1):(3 * p)], r_oos %*% w_bootparam[i, (2*p + 1):(3 * p)], r_oos %*% w_boot[i, (2*p + 1):(3 * p)], r_oos %*% w_factor_bootparam_pca[i, (2*p + 1):(3 * p)], r_oos %*% w_factor_boot_pca[i, (2*p + 1):(3 * p)], r_oos %*% w_factor_bootparam_obs[i, (2*p + 1):(3 * p)], r_oos %*% w_factor_boot_obs[i, (2*p + 1):(3 * p)])
}



# Saving Results
colnames(Rport) <- c("mv_Markowitz", "mv_MichaudParam", "mv_MichaudNonP", "mv_FactorParamPCA", "mv_FactorNonPPCA", "mv_FactorParamObs", "mv_FactorNonPObs",
                    "tp_Markowitz", "tp_MichaudParam", "tp_MichaudNonP", "tp_FactorParamPCA", "tp_FactorNonPPCA", "tp_FactorParamObs", "tp_FactorNonPObc",
                    "ef_Markowitz", "ef_MichaudParam", "ef_MichaudNonP", "ef_FactorParamPCA", "ef_FactorNonPPCA", "ef_FactorParamObs", "ef_FactorNonPObc")
write.table(Rport, paste0("Results/Rport_", InS, "_ibov.csv"), sep = ",")


colnames(to) <- c("mv_Markowitz", "tp_Markowitz", "ef_Markowitz",
                     "mv_MichaudParam", "tp_MichaudParam", "ef_MichaudParam",
                     "mv_MichaudNonP", "tp_MichaudNonP", "ef_MichaudNonP",
                     "mv_FactorParamPCA", "tp_FactorParamPCA", "ef_FactorParamPCA",
                     "mv_FactorNonPPCA", "tp_FactorNonPPCA", "ef_FactorNonPPCA",
                     "mv_FactorParamFF", "tp_FactorParamFF", "ef_FactorParamFF",
                     "mv_FactorNonPFF", "tp_FactorNonPFF",  "ef_FactorNonPFF")
write.table(to, paste0("Results/to_", InS, "_ibov.csv"), sep = ",")


colnames(sspw) <- c("mv_Markowitz", "tp_Markowitz", "ef_Markowitz",
                    "mv_MichaudParam", "tp_MichaudParam", "ef_MichaudParam",
                    "mv_MichaudNonP", "tp_MichaudNonP", "ef_MichaudNonP",
                    "mv_FactorParamPCA", "tp_FactorParamPCA", "ef_FactorParamPCA",
                    "mv_FactorNonPPCA", "tp_FactorNonPPCA", "ef_FactorNonPPCA",
                    "mv_FactorParamFF", "tp_FactorParamFF", "ef_FactorParamFF",
                    "mv_FactorNonPFF", "tp_FactorNonPFF",  "ef_FactorNonPFF")
write.table(sspw, paste0("Results/sspw_", InS, "_ibov.csv"), sep = ",")



# Tables for short-selling MVP
Caption <- "Out-of-sample performance measures of the minimum variance portfolio with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,1:nmethods], 2, medidas, rf), 
                     apply(to[, seq(1, nconstrains*nmethods, by = 3)], 2, mean), 
                     apply(sspw[, seq(1, nconstrains*nmethods, by = 3)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_mvp") %>% 
  save_kable(keep_tex = T, file = paste0("Results/mvp_", InS, "_ibov.tex"))



# Tables for short-selling TP
Caption <- "Out-of-sample performance measures of the tangency portfolio with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,(nmethods + 1):(2*nmethods)], 2, medidas, rf),  
                     apply(to[, seq(2, nconstrains*nmethods, by = 3)], 2, mean), 
                     apply(sspw[, seq(2, nconstrains*nmethods, by = 3)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_tp") %>% 
  save_kable(keep_tex = T, file = paste0("Results/tp_", InS, "_ibov.tex"))

# Tables for short-selling EF
Caption <- "Out-of-sample performance measures of the eficient frontier with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,(2*nmethods + 1):(3*nmethods)], 2, medidas, rf),  
                     apply(to[, seq(3, nconstrains*nmethods, by = 3)], 2, mean), 
                     apply(sspw[, seq(3, nconstrains*nmethods, by = 3)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_ef") %>% 
  save_kable(keep_tex = T, file = paste0("Results/ef_", InS, "_ibov.tex"))

