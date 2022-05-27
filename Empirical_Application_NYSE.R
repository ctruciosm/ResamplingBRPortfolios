#################################
###   Empirical Application   ###
#################################
rm(list = ls())
library(dplyr)
library(mvtnorm)
library(moments)
library(RiskPortfolios)
library(fPortfolio)
library(tidyr)
library(readxl)
library(stringr)
library(dplyr)
library(lubridate)
library(POET)
library(Matrix)
library(kableExtra)
library(nortsTest)
source("Resampling_Techniques.R")
source("Performance_Measures.R")
source("Auxiliary_Functions.R")

# Importing and Wrangling  Data
{
  monthly_data <- read_xlsx("Data/economatica_nyse.xlsx",  skip = 3)
  colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\ndo fechamento\nem 1 mês\nEm moeda orig\najust p/ prov\n", "")
  
  monthly_data <- monthly_data %>%
    mutate(Data = str_replace(Data, "Jan", "01"), 
           Data = str_replace(Data, "Fev", "02"),
           Data = str_replace(Data, "Mar", "03"), 
           Data = str_replace(Data, "Abr", "04"),
           Data = str_replace(Data, "Mai", "05"), 
           Data = str_replace(Data, "Jun", "06"),
           Data = str_replace(Data, "Jul", "07"), 
           Data = str_replace(Data, "Ago", "08"),
           Data = str_replace(Data, "Set", "09"), 
           Data = str_replace(Data, "Out", "10"),
           Data = str_replace(Data, "Nov", "11"), 
           Data = str_replace(Data, "Dez", "12")) %>% 
    mutate(Data = lubridate::my(Data)) %>% 
    mutate_if(is.character, as.numeric) %>% 
    filter(Data >= '1995-06-01', Data <= '2022-01-01')
  
  monthly_data_dates <- monthly_data %>% 
    select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) 
  
  monthly_data <- monthly_data_dates %>% select(-Data)
  
  famafrench <- read.csv("Data/F_F_Research_Data_Factors.CSV") %>% 
    mutate(Date = lubridate::ym(Date)) %>% 
    filter(Date >= '1995-06-01', Date <= '2022-01-01') %>% 
    select(-Date)
  famafrench5f <- read.csv("Data/F_F_Research_Data_5_Factors_2x3.csv") %>% 
    mutate(Date = lubridate::ym(Date)) %>% 
    filter(Date >= '1995-06-01', Date <= '2022-01-01') %>% 
    select(-Date)
}

monthly_data_longer <- pivot_longer(monthly_data, cols = everything(), values_to = "returns", names_to = "assets")


assets_names <- 
  monthly_data_longer %>% 
  group_by(assets) %>% 
  summarise(LjungBox = Box.test(returns, type =  "Ljung-Box")$p.value,
            LjungBox2 = Box.test(returns^2, type =  "Ljung-Box")$p.value) %>% 
  filter(LjungBox > 0.05, LjungBox2 > 0.05) %>% 
  select(assets)

monthly_data <-  monthly_data %>% select(as.vector(as.matrix(assets_names)))

InS <- 120
OoS <- nrow(monthly_data) - InS
p <- ncol(monthly_data) 
nboot <- 1000
nmethods <- 9

#######################################
###   Minimum Variance Portfolios   ###
#######################################
nconstrains <- 3
w_estim = w_bootparam = w_boot = w_factor_bootparam_pca = w_factor_boot_pca = w_factor_bootparam_ff = w_factor_boot_ff = w_factor_bootparam_ff5f = w_factor_boot_ff5f = matrix(NA, nrow = OoS, ncol = nconstrains*p)
Rport <- matrix(NA, nrow = OoS, ncol = nconstrains*nmethods)
w_measures <- matrix(0, nrow = OoS, ncol = 2*nconstrains*nmethods)
for (i in 1:OoS) {
  print(i)
  
  r_ins <- as.matrix(monthly_data[i:(InS - 1 + i), ])
  ff <- as.matrix(famafrench[i:(InS - 1 + i), ])
  ff5f <- as.matrix(famafrench5f[i:(InS - 1 + i), ])
  observed_factors <- list(ff, ff5f)
  r_oos <- as.numeric(monthly_data[InS - 1 + i + 1, ]) 
  
  set.seed(i + 1234)
  weights_unc <- calculate_portfolio_weights(x = r_ins, constrains_opt = list(type = 'minvol'), nboot, factors = observed_factors)
  weights_ssc <- calculate_portfolio_weights(x = r_ins, constrains_opt = list(type = 'minvol', constraint = 'lo'), nboot, factors = observed_factors)
  weights_luc <- calculate_portfolio_weights(x = r_ins, constrains_opt = list(type = 'minvol', constraint = 'user', LB = rep(0, p), UB = rep(3/p, p)), nboot, factors = observed_factors)
  
  w_estim[i, ] <- c(weights_unc[1, ], weights_ssc[1, ], weights_luc[1, ])
  w_bootparam[i, ] <- c(weights_unc[2, ], weights_ssc[2, ], weights_luc[2, ])
  w_boot[i, ] <- c(weights_unc[3, ], weights_ssc[3, ], weights_luc[3, ])
  w_factor_bootparam_pca[i, ] <- c(weights_unc[4, ], weights_ssc[4, ], weights_luc[4, ])
  w_factor_boot_pca[i, ] <- c(weights_unc[5, ], weights_ssc[5, ], weights_luc[5, ])
  w_factor_bootparam_ff[i, ] <- c(weights_unc[6, ], weights_ssc[6, ], weights_luc[6, ])
  w_factor_boot_ff[i, ] <- c(weights_unc[7, ], weights_ssc[7, ], weights_luc[7, ])
  w_factor_bootparam_ff5f[i, ] <- c(weights_unc[8, ], weights_ssc[8, ], weights_luc[8, ])
  w_factor_boot_ff5f[i, ] <- c(weights_unc[9, ], weights_ssc[9, ], weights_luc[9, ])
  
   
  if (i > 1) {
    w_measures[i, ] <- c(unlist(weights_measures(w_estim[i - 1, 1:p], w_estim[i, 1:p], r_oos)),
                         unlist(weights_measures(w_bootparam[i - 1, 1:p], w_bootparam[i, 1:p], r_oos)),
                         unlist(weights_measures(w_boot[i - 1, 1:p], w_boot[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_pca[i - 1, 1:p], w_factor_bootparam_pca[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_boot_pca[i - 1, 1:p], w_factor_boot_pca[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff[i - 1, 1:p], w_factor_bootparam_ff[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff[i - 1, 1:p], w_factor_boot_ff[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff5f[i - 1, 1:p], w_factor_bootparam_ff5f[i, 1:p], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff5f[i - 1, 1:p], w_factor_boot_ff5f[i, 1:p], r_oos)),
                         unlist(weights_measures(w_estim[i - 1, (p + 1):(2*p)], w_estim[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_bootparam[i - 1, (p + 1):(2*p)], w_bootparam[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_boot[i - 1, (p + 1):(2*p)], w_boot[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_pca[i - 1, (p + 1):(2*p)], w_factor_bootparam_pca[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_pca[i - 1, (p + 1):(2*p)], w_factor_boot_pca[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff[i - 1, (p + 1):(2*p)], w_factor_bootparam_ff[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff[i - 1, (p + 1):(2*p)], w_factor_boot_ff[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff5f[i - 1, (p + 1):(2*p)], w_factor_bootparam_ff5f[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff5f[i - 1, (p + 1):(2*p)], w_factor_boot_ff5f[i, (p + 1):(2*p)], r_oos)),
                         unlist(weights_measures(w_estim[i - 1, (2*p + 1):(3*p)], w_estim[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_bootparam[i - 1, (2*p + 1):(3*p)], w_bootparam[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_boot[i - 1, (2*p + 1):(3*p)], w_boot[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_pca[i - 1, (2*p + 1):(3*p)], w_factor_bootparam_pca[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_pca[i - 1, (2*p + 1):(3*p)], w_factor_boot_pca[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff[i - 1, (2*p + 1):(3*p)], w_factor_bootparam_ff[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff[i - 1, (2*p + 1):(3*p)], w_factor_boot_ff[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_bootparam_ff5f[i - 1, (2*p + 1):(3*p)], w_factor_bootparam_ff5f[i, (2*p + 1):(3*p)], r_oos)),
                         unlist(weights_measures(w_factor_boot_ff5f[i - 1, (2*p + 1):(3*p)], w_factor_boot_ff5f[i, (2*p + 1):(3*p)], r_oos)))
  }
  Rport[i,] <-  c(r_oos %*% w_estim[i, 1:p], r_oos %*% w_bootparam[i, 1:p], r_oos %*% w_boot[i, 1:p], r_oos %*% w_factor_bootparam_pca[i, 1:p], r_oos %*% w_factor_boot_pca[i, 1:p], r_oos %*% w_factor_bootparam_ff[i, 1:p], r_oos %*% w_factor_boot_ff[i, 1:p], r_oos %*% w_factor_bootparam_ff5f[i, 1:p], r_oos %*% w_factor_boot_ff5f[i, 1:p],
                 r_oos %*% w_estim[i, (p + 1):(2 * p)], r_oos %*% w_bootparam[i, (p + 1):(2 * p)], r_oos %*% w_boot[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam_pca[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot_pca[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam_ff[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot_ff[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam_ff5f[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot_ff5f[i, (p + 1):(2 * p)],
                 r_oos %*% w_estim[i, (2 * p + 1):(3 * p)], r_oos %*% w_bootparam[i, (2 * p + 1):(3 * p)], r_oos %*% w_boot[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_bootparam_pca[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_boot_pca[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_bootparam_ff[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_boot_ff[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_bootparam_ff5f[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_boot_ff5f[i, (2 * p + 1):(3 * p)])
}

# Saving Results
colnames(Rport) <- c("unc_Markowitz", "unc_MichaudParam", "unc_MichaudNonP", "unc_FactorParamPCA", "unc_FactorNonPPCA", "unc_FactorParamFF", "unc_FactorNonPFF", "unc_FactorParamFF5f", "unc_FactorNonPFF5f",
                     "ssc_Markowitz", "ssc_MichaudParam", "ssc_MichaudNonP", "ssc_FactorParamPCA", "ssc_FactorNonPPCA", "ssc_FactorParamFF", "ssc_FactorNonPFF", "ssc_FactorParamFF5f", "ssc_FactorNonPFF5f",
                     "luc_Markowitz", "luc_MichaudParam", "luc_MichaudNonP", "luc_FactorParamPCA", "luc_FactorNonPPCA", "luc_FactorParamFF", "luc_FactorNonPFF", "luc_FactorParamFF5f", "luc_FactorNonPFF5f")
write.table(Rport, paste0("Results/Rport_", InS, "_nyse.csv"), sep = ",")


# Tables for unscontrained MVP
Caption <- "Out-of-sample performance measures of the unconstrained minimum variance portfolio: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,1:nmethods],2,medidas), apply(w_measures[-1, seq(1, 2*nmethods, by = 2)], 2, mean),apply(w_measures[-1, seq(2, 2*nmethods, by = 2)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Fama-French", "Factor-Based Non-Parametric Fama-French", "Factor-Based Parametric Fama-French 5F", "Factor-Based Non-Parametric Fama-French 5F")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "unc_mvp") %>% 
  save_kable(keep_tex = T, file = paste0("Results/unc_mvp_", InS, "_nyse.tex"))



# Tables for short-selling MVP
Caption <- "Out-of-sample performance measures of the minimum variance portfolio with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,(nmethods + 1):(2*nmethods)],2,medidas), apply(w_measures[-1, seq(2*nmethods + 1, 4*nmethods, by = 2)], 2, mean),apply(w_measures[-1, seq(2*nmethods + 2, 4*nmethods, by = 2)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Fama-French", "Factor-Based Non-Parametric Fama-French", "Factor-Based Parametric Fama-French 5F", "Factor-Based Non-Parametric Fama-French 5F")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_mvp") %>% 
  save_kable(keep_tex = T, file = paste0("Results/ssc_mvp_", InS, "_nyse.tex"))

# Tables for lower-upper MVP
Caption <- "Out-of-sample performance measures of the minimum variance portfolio with lower (0) and upper (10) bound constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,(2*nmethods + 1):(3*nmethods)],2,medidas), apply(w_measures[-1, seq(4*nmethods + 1, 6*nmethods, by = 2)], 2, mean),apply(w_measures[-1, seq(4*nmethods + 2, 6*nmethods, by = 2)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric","Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Factor-Based Parametric Fama-French", "Factor-Based Non-Parametric Fama-French", "Factor-Based Parametric Fama-French 5F", "Factor-Based Non-Parametric Fama-French 5F")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "luc_mvp") %>% 
  save_kable(keep_tex = T, file = paste0("Results/luc_mvp_", InS, "_nyse.tex"))

