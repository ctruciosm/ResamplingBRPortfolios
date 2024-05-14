#################################
###   Empirical Application   ###
#################################
#setwd("/home/ctrucios/Dropbox/Research/Resampling")
rm(list = ls())
library(mvtnorm)
library(POET)
library(Matrix)
library(kableExtra)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(lubridate)
library(nortsTest)
source("Resampling_Techniques.R")
source("Auxiliary_Functions.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#            Data              #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
monthly_data <- read_xlsx("economatica_nyse_mensal.xlsx", na = c("-", "NA"), skip = 3, col_types = c("text", rep("numeric", 1419)))
colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\ndo fechamento\nem 1 meses\nEm moeda orig\najust p/ prov\n", "")

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
  dplyr::filter(Data >= '1995-01-01', Data < '2023-01-01')
monthly_data <- monthly_data %>% dplyr::select(-Data)

ff <- read.csv("F-F_Research_Data_Factors.csv") |> 
  mutate(Data = lubridate::ym(Data)) |> 
  dplyr::filter(Data >= '1995-01-01', Data < '2023-01-01') |> 
  select(-Data)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#       General Settings       #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
InS <- 120
OoS <- nrow(monthly_data) - InS
p <- ncol(monthly_data) 
nboot <- 2000
nmethods <- 10
rf <- 0.5
lambda <- 2

Rport <- matrix(NA, nrow = OoS, ncol = 2*nmethods)
to <- matrix(0, nrow = OoS - 1, ncol = 2*nmethods)
sspw <- matrix(0, nrow = OoS, ncol = 2*nmethods)

w_ew = matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(monthly_data)))
w_mv_estim = w_mv_bootparam = w_mv_boot = w_mv_estim_pca = w_mv_factor_bootparam_pca = w_mv_factor_boot_pca = w_mv_estim_obs = w_mv_factor_bootparam_obs = w_mv_factor_boot_obs = matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(monthly_data)))
w_ef_estim = w_ef_bootparam = w_ef_boot = w_ef_estim_pca = w_ef_factor_bootparam_pca = w_ef_factor_boot_pca = w_ef_estim_obs = w_ef_factor_bootparam_obs = w_ef_factor_boot_obs = matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(monthly_data)))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#    Out-of-sample exercise    #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
for (i in 1:OoS) { 
  print(i)
  
  r_aux <- monthly_data[i:(InS + i), ] %>% select(names(which(apply(is.na(monthly_data[i:(InS + i), ]), 2, sum) == 0)))
  r_ins <- as.matrix(r_aux[1:InS, ])
  ff_ins <- as.matrix(ff[i:(InS - 1 + i), ])
  r_oos <- as.numeric(monthly_data[InS + i, colnames(r_ins)]) 
  names_window <- names(r_aux)
  
  set.seed(i + 2468)
  w_mv_aux <- calculate_portfolio_weights(x = r_ins, type = "mv",  nboot, factors = ff_ins)
  w_ef_aux <- calculate_portfolio_weights(x = r_ins, type = "ef",  nboot, factors = ff_ins, lambda = lambda)
  
  col_window <- ncol(r_ins)
  w_ew[i, names_window] = rep(1/col_window, col_window)
  
  w_mv_estim[i, names_window] <- w_mv_aux$weights[1, ]
  w_ef_estim[i, names_window] <- w_ef_aux$weights[1, ]
  
  w_mv_bootparam[i, names_window] <- w_mv_aux$weights[2, ]
  w_ef_bootparam[i, names_window] <- w_ef_aux$weights[2, ]
  
  w_mv_boot[i, names_window] <- w_mv_aux$weights[3, ]
  w_ef_boot[i, names_window] <- w_ef_aux$weights[3, ]
  
  w_mv_estim_pca[i, names_window] <- w_mv_aux$weights[4, ]
  w_ef_estim_pca[i, names_window] <- w_ef_aux$weights[4, ]
  
  w_mv_factor_bootparam_pca[i, names_window] <- w_mv_aux$weights[5, ]
  w_ef_factor_bootparam_pca[i, names_window] <- w_ef_aux$weights[5, ]
  
  w_mv_factor_boot_pca[i, names_window] <- w_mv_aux$weights[6, ]
  w_ef_factor_boot_pca[i, names_window] <- w_ef_aux$weights[6, ]
  
  w_mv_estim_obs[i, names_window] <- w_mv_aux$weights[7, ]
  w_ef_estim_obs[i, names_window] <- w_ef_aux$weights[7, ]
  
  w_mv_factor_bootparam_obs[i, names_window] <- w_mv_aux$weights[8, ]
  w_ef_factor_bootparam_obs[i, names_window] <- w_ef_aux$weights[8, ]
  
  w_mv_factor_boot_obs[i, names_window] <- w_mv_aux$weights[9, ]
  w_ef_factor_boot_obs[i, names_window] <- w_ef_aux$weights[9, ]
  
  sspw[i, ] <- c(sum(w_ew[i, ]^2), sum(w_mv_estim[i, names_window]^2), sum(w_mv_bootparam[i, names_window]^2), sum(w_mv_boot[i, names_window]^2), sum(w_mv_estim_pca[i, names_window]^2), sum(w_mv_factor_bootparam_pca[i, names_window]^2), sum(w_mv_estim_obs[i, names_window]^2),  sum(w_mv_factor_boot_pca[i, names_window]^2), sum(w_mv_factor_bootparam_obs[i, names_window]^2), sum(w_mv_factor_boot_obs[i, names_window]^2),
                 sum(w_ew[i, ]^2), sum(w_ef_estim[i, names_window]^2), sum(w_ef_bootparam[i, names_window]^2), sum(w_ef_boot[i, names_window]^2), sum(w_ef_estim_pca[i, names_window]^2), sum(w_ef_factor_bootparam_pca[i, names_window]^2), sum(w_ef_estim_obs[i, names_window]^2),  sum(w_ef_factor_boot_pca[i, names_window]^2), sum(w_ef_factor_bootparam_obs[i, names_window]^2), sum(w_ef_factor_boot_obs[i, names_window]^2))
  
  if (i > 1) {
    to[i - 1, ] <- c(calculate_turnover(w_ew[i - 1,], w_ew[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_estim[i - 1,], w_mv_estim[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_bootparam[i - 1,], w_mv_bootparam[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_boot[i - 1,], w_mv_boot[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_estim_pca[i - 1,], w_mv_estim_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_factor_bootparam_pca[i - 1,], w_mv_factor_bootparam_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_factor_boot_pca[i - 1,], w_mv_factor_boot_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_estim_obs[i - 1,], w_mv_estim_obs[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_factor_bootparam_obs[i - 1,], w_mv_factor_bootparam_obs[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_mv_factor_boot_obs[i - 1,], w_mv_factor_boot_obs[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ew[i - 1,], w_ew[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_estim[i - 1,], w_ef_estim[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_bootparam[i - 1,], w_ef_bootparam[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_boot[i - 1,], w_ef_boot[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_estim_pca[i - 1,], w_ef_estim_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_factor_bootparam_pca[i - 1,], w_ef_factor_bootparam_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_factor_boot_pca[i - 1,], w_ef_factor_boot_pca[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_estim_obs[i - 1,], w_ef_estim_obs[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_factor_bootparam_obs[i - 1,], w_ef_factor_bootparam_obs[i,], monthly_data[InS + i, ]),
                     calculate_turnover(w_ef_factor_boot_obs[i - 1,], w_ef_factor_boot_obs[i,], monthly_data[InS + i, ]))
  }
  
  Rport[i,] <-  c(mean(r_oos), r_oos %*% w_mv_estim[i, names_window], r_oos %*% w_mv_bootparam[i, names_window], r_oos %*% w_mv_boot[i, names_window],  r_oos %*% w_mv_estim_pca[i, names_window],  r_oos %*% w_mv_factor_bootparam_pca[i, names_window], r_oos %*% w_mv_factor_boot_pca[i, names_window], r_oos %*% w_mv_estim_obs[i, names_window],  r_oos %*% w_mv_factor_bootparam_obs[i, names_window], r_oos %*% w_mv_factor_boot_obs[i, names_window],
                  mean(r_oos), r_oos %*% w_ef_estim[i, names_window], r_oos %*% w_ef_bootparam[i, names_window], r_oos %*% w_ef_boot[i, names_window], r_oos %*% w_ef_estim_pca[i, names_window],  r_oos %*% w_ef_factor_bootparam_pca[i, names_window], r_oos %*% w_ef_factor_boot_pca[i, names_window], r_oos %*% w_ef_estim_obs[i, names_window], r_oos %*% w_ef_factor_bootparam_obs[i, names_window],  r_oos %*% w_ef_factor_boot_obs[i, names_window])
  
}

# Saving Results

colnames(Rport) <- c("ew",
                     "mv_Markowitz", "mv_MichaudParam", "mv_MichaudNonP", 
                     "mv_MarkowitzPCA","mv_FactorParamPCA", "mv_FactorNonPPCA", 
                     "mv_MarkowitzObs","mv_FactorParamObs", "mv_FactorNonPObs",
                     "ef_ew",
                     "ef_Markowitz", "ef_MichaudParam", "ef_MichaudNonP", 
                     "ef_MarkowitzPCA", "ef_FactorParamPCA", "ef_FactorNonPPCA", 
                     "ef_MarkowitzObs", "ef_FactorParamObs", "ef_FactorNonPObc")
write.table(Rport, paste0("Results/Rport_2_", InS, "_nyse.csv"), sep = ",")


colnames(to) <- c("ew",
                  "mv_Markowitz", "mv_MichaudParam", "mv_MichaudNonP", 
                  "mv_MarkowitzPCA","mv_FactorParamPCA", "mv_FactorNonPPCA", 
                  "mv_MarkowitzObs","mv_FactorParamObs", "mv_FactorNonPObs",
                  "ef_ew",
                  "ef_Markowitz", "ef_MichaudParam", "ef_MichaudNonP", 
                  "ef_MarkowitzPCA", "ef_FactorParamPCA", "ef_FactorNonPPCA", 
                  "ef_MarkowitzObs", "ef_FactorParamObs", "ef_FactorNonPObc")
write.table(to, paste0("Results/to_2_", InS, "_nyse.csv"), sep = ",")


colnames(sspw) <- c("ew",
                    "mv_Markowitz", "mv_MichaudParam", "mv_MichaudNonP", 
                    "mv_MarkowitzPCA","mv_FactorParamPCA", "mv_FactorNonPPCA", 
                    "mv_MarkowitzObs","mv_FactorParamObs", "mv_FactorNonPObs",
                    "ef_ew",
                    "ef_Markowitz", "ef_MichaudParam", "ef_MichaudNonP", 
                    "ef_MarkowitzPCA", "ef_FactorParamPCA", "ef_FactorNonPPCA", 
                    "ef_MarkowitzObs", "ef_FactorParamObs", "ef_FactorNonPObc")
write.table(sspw, paste0("Results/sspw_2_", InS, "_nyse.csv"), sep = ",")



# Tables for short-selling MVP
Caption <- "Out-of-sample performance measures of the minimum variance portfolio with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport, 2, medidas, rf), apply(to, 2, mean), apply(sspw, 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("EW", "Markowitz", "Michaud Parametric", "Michaud Non-Parametric", "Markowitz-PCA", "Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Markowitz-Ibov", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>%
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_mvp") %>%
  save_kable(keep_tex = T, file = paste0("Results/mvp_2_", InS, "_nyse.tex"))



# Tables for short-selling EF
Caption <- "Out-of-sample performance measures of the eficient frontier with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport, 2, medidas, rf), apply(to, 2, mean), apply(sspw, 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("EW","Markowitz", "Michaud Parametric", "Michaud Non-Parametric", "Markowitz-PCA", "Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Markowitz-Ibov", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "ssc_ef") %>% 
  save_kable(keep_tex = T, file = paste0("Results/ef_2_", InS, "_nyse.tex"))

