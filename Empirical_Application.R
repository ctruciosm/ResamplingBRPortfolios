#################################
###   Empirical Application   ###
#################################
rm(list = ls())
library(dplyr)
library(mvtnorm)
library(RiskPortfolios)
library(tidyr)
library(readxl)
library(stringr)
library(dplyr)
library(lubridate)
library(POET)
library(Matrix)
source("Resampling_Techniques.R")
source("Performance_Measures.R")
source("Auxiliary_Functions.R")

# Importing and Wrangling  Data
{
monthly_data <- read_xlsx("Data/dados_mensais_IBRX.xlsx", na = "-")
colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\r\ndo fechamento\r\nem 1 mês\r\nEm moeda orig\r\najust p/ prov\r\n", "")

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
  filter(Data >= '2000-01-01')
monthly_data <- monthly_data %>% select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) 
monthly_data <- monthly_data[, -1]
}

InS <- 60
OoS <- nrow(monthly_data) - InS
p <- ncol(monthly_data) 
nboot <- 250
nmethods <- 5

#######################################
###   Minimum Variance Portfolios   ###
#######################################
nconstrains <- 3
w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = matrix(NA, nrow = OoS, ncol = nconstrains*p)
Rport <- matrix(NA, nrow = OoS, ncol = nconstrains*nmethods)
for (i in 1:OoS) {
  print(i)
  set.seed(i + 531)
  r_ins <- as.matrix(monthly_data[i:(InS - 1 + i), ])
  r_oos <- as.numeric(monthly_data[InS - 1 + i + 1, ]) 
  
  weights_unc <- calculate_portfolio_weights(r_ins, constrains_opt = list(type = 'minvol'), nboot)
  weights_ssc <- calculate_portfolio_weights(r_ins, constrains_opt = list(type = 'minvol', constraint = 'lo'), nboot)
  weights_luc <- calculate_portfolio_weights(r_ins, constrains_opt = list(type = 'minvol', constraint = 'user', LB = rep(0, p), UB = rep(0.1, p)), nboot)
  
  w_estim[i, ] <- c(weights_unc[1, ], weights_ssc[1, ], weights_luc[1, ])
  w_bootparam[i, ] <- c(weights_unc[2, ], weights_ssc[2, ], weights_luc[2, ])
  w_boot[i, ] <- c(weights_unc[3, ], weights_ssc[3, ], weights_luc[3, ])
  w_factor_bootparam[i, ] <- c(weights_unc[4, ], weights_ssc[4, ], weights_luc[4, ])
  w_factor_boot[i, ] <- c(weights_unc[5, ], weights_ssc[5, ], weights_luc[5, ])

  Rport[i,] <-  c(r_oos %*% w_estim[i, 1:p], r_oos %*% w_bootparam[i, 1:p], r_oos %*% w_boot[i, 1:p], r_oos %*% w_factor_bootparam[i, 1:p], r_oos %*% w_factor_boot[i, 1:p],
                 r_oos %*% w_estim[i, (p + 1):(2 * p)], r_oos %*% w_bootparam[i, (p + 1):(2 * p)], r_oos %*% w_boot[i, (p + 1):(2 * p)], r_oos %*% w_factor_bootparam[i, (p + 1):(2 * p)], r_oos %*% w_factor_boot[i, (p + 1):(2 * p)],
                 r_oos %*% w_estim[i, (2 * p + 1):(3 * p)], r_oos %*% w_bootparam[i, (2 * p + 1):(3 * p)], r_oos %*% w_boot[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_bootparam[i, (2 * p + 1):(3 * p)], r_oos %*% w_factor_boot[i, (2 * p + 1):(3 * p)])

}

oos_results <- apply(Rport[,1:nmethods],2,medidas)
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric",
                            "Factor-Based Parametric", "Factor-Based Non-Parametric")

xtable::xtable(t(oos_results),digits = 4)




