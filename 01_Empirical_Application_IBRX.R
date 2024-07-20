#################################
###   Empirical Application   ###
#################################
rm(list = ls())
library(mvtnorm)
library(moments)
library(POET)
library(Matrix)
library(quadprog)
library(readxl)
library(kableExtra)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
source("Resampling_Techniques.R")
source("Auxiliary_Functions.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#            Data              #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
composition <- read_excel("Data/composicao_IBRx.xlsx") |> select(-Company, -Type) |> 
  pivot_longer(cols = `Dec-97`:`May-24`, names_to = "Date", values_to = "valores") |> 
  mutate(dates = lubridate::my(Date)) |>
  select(-Date) |> 
  select(dates, everything()) |> 
  pivot_wider(names_from = Code, values_from = valores)

stocks <- read_excel("Data/economatica_b3.xlsx", skip = 3, na = c("-", "NA")) |> 
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
         Data = str_replace(Data, "Dez", "12")) |> 
  mutate(dates = lubridate::my(Data)) |> 
  select(-Data) |> 
  select(dates, everything()) |> 
  filter(dates >= '1997-12-01') |> 
  filter(dates <= '2023-12-01')
colnames(stocks) <- str_replace(colnames(stocks), "Retorno\ndo fechamento\nem 1 meses\nEm moeda orig\najust p/ prov\n", "")

ibovespa <- read_excel("Data/economatica_ibov.xlsx", skip = 3, na = c("-", "NA")) |> 
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
         Data = str_replace(Data, "Dez", "12")) |> 
  mutate(dates = lubridate::my(Data)) |> 
  select(-Data) |> 
  select(dates, everything()) |> 
  filter(dates >= '1997-12-01') |> 
  filter(dates <= '2023-12-01')
colnames(ibovespa) <- str_replace(colnames(ibovespa), "Retorno\ndo fechamento\nem 1 meses\nEm moeda orig\najust p/ prov\n", "")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#       General Settings       #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
InS <- 120
OoS <- nrow(stocks) - InS
p <- ncol(stocks) - 1
nboot <- 500
rf <- 0.5
nmethods <- 10
lambda <- 5

w_estim_mv = w_bootparam_mv = w_boot_mv = w_estim_pca_mv = w_factor_bootparam_pca_mv = w_factor_boot_pca_mv = w_estim_obs_mv = w_factor_bootparam_obs_mv = w_factor_boot_obs_mv = matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(stocks)[-1]))
w_estim_ef = w_bootparam_ef = w_boot_ef = w_estim_pca_ef = w_factor_bootparam_pca_ef = w_factor_boot_pca_ef = w_estim_obs_ef = w_factor_bootparam_obs_ef = w_factor_boot_obs_ef = matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(stocks)[-1]))
w_ew_full <- matrix(0, nrow = OoS, ncol = p, dimnames = list(NULL, colnames(stocks)[-1]))

Rport <- matrix(NA, nrow = OoS, ncol = 2*nmethods)
to <- matrix(0, nrow = OoS - 1, ncol = 2*nmethods)
sspw <- matrix(0, nrow = OoS, ncol = 2*nmethods)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#    Out-of-sample exercise    #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
for (i in 1:OoS) { 
  print(i)
  date_ins <- stocks[i:(InS - 1 + i), ]$dates
  aux_1 <- composition[i:(InS - 1 + i), ] |> filter(dates <= date_ins[InS]) |> tail(1) |> select(-dates) |> is.na() |> apply(2, sum) 
  aux_2 <- stocks[i:(InS - 1 + i), ] |> select(-dates) |> is.na() |> apply(2, sum) 
  aux <- intersect(names(which(aux_1 == 0)), names(which(aux_2 == 0)))

  retu_ins <- stocks[i:(InS - 1 + i), ] |> select(aux) |> as.matrix()
  r_oos <- stocks[InS + i, ] |> select(aux) |> as.matrix()
  r_oos_full <- matrix(NA, ncol = p, nrow = 1, dimnames = list(NULL, colnames(stocks)[-1]))
  r_oos_full[1, aux] <- r_oos
  ibov_ins <- ibovespa[i:(InS - 1 + i), ] |> select(IBOV) |> as.matrix()
  w_ew = rep(1/length(aux), length(aux))
    
  set.seed(i)
  #w_mv <- calculate_portfolio_weights(x = retu_ins, type = "mv",  nboot, factors = ibov_ins)
  w_ef <- calculate_portfolio_weights(x = retu_ins, type = "ef",  nboot, factors = ibov_ins, lambda = lambda)
  w_mv <- w_ef
  
  Rport[i,] <- c(mean(r_oos), w_mv %*% t(r_oos), mean(r_oos), w_ef %*% t(r_oos))
  sspw[i, ] <- as.numeric(c(sum(w_ew^2), apply(w_mv^2, 1, sum), sum(w_ew^2), apply(w_ef^2, 1, sum)))
    
  w_ew_full[i, aux] <- rep(1/length(aux), length(aux))
  w_estim_mv[i, aux] <- w_mv[1, ]
  w_estim_ef[i, aux] <- w_ef[1, ]
  w_bootparam_mv[i, aux] <- w_mv[2, ]
  w_bootparam_ef[i, aux] <- w_ef[2, ]
  w_boot_mv[i, aux] <- w_mv[3, ]
  w_boot_ef[i, aux] <- w_ef[3, ]
  w_estim_pca_mv[i, aux] <- w_mv[4, ]
  w_estim_pca_ef[i, aux] <- w_ef[4, ]
  w_factor_bootparam_pca_mv[i, aux] <- w_mv[5, ]
  w_factor_bootparam_pca_ef[i, aux] <- w_ef[5, ]
  w_factor_boot_pca_mv[i, aux] <- w_mv[6, ]
  w_factor_boot_pca_ef[i, aux] <- w_ef[6, ]
  w_estim_obs_mv[i, aux] <- w_mv[7, ]
  w_estim_obs_ef[i, aux] <- w_ef[7, ]
  w_factor_bootparam_obs_mv[i, aux] <- w_mv[8, ]
  w_factor_bootparam_obs_ef[i, aux] <- w_ef[8, ]
  w_factor_boot_obs_mv[i,aux] <- w_mv[9, ]
  w_factor_boot_obs_ef[i,aux] <- w_ef[9, ]
   
  if (i > 2) {
    to[i - 1, ] <- c(calculate_to(w_ew_full[i - 1, ], w_ew_full[i, ], r_oos_full, p),
                     calculate_to(w_estim_mv[i - 1, ], w_estim_mv[i, ], r_oos_full, p),
                     calculate_to(w_bootparam_mv[i - 1, ], w_bootparam_mv[i, ], r_oos_full, p),
                     calculate_to(w_boot_mv[i - 1, ], w_boot_mv[i, ], r_oos_full, p),
                     calculate_to(w_estim_pca_mv[i - 1, ], w_estim_pca_mv[i, ], r_oos_full, p),
                     calculate_to(w_factor_bootparam_pca_mv[i - 1, ], w_factor_bootparam_pca_mv[i, ], r_oos_full, p),
                     calculate_to(w_factor_boot_pca_mv[i - 1, ], w_factor_boot_pca_mv[i, ], r_oos_full, p),
                     calculate_to(w_estim_obs_mv[i - 1, ], w_estim_obs_mv[i, ], r_oos_full, p),
                     calculate_to(w_factor_bootparam_obs_mv[i - 1, ], w_factor_bootparam_obs_mv[i, ], r_oos_full, p),
                     calculate_to(w_factor_boot_obs_mv[i - 1, ], w_factor_boot_obs_mv[i, ], r_oos_full, p),
                     calculate_to(w_ew_full[i - 1, ], w_ew_full[i, ], r_oos_full, p),
                     calculate_to(w_estim_ef[i - 1, ], w_estim_ef[i, ], r_oos_full, p),
                     calculate_to(w_bootparam_ef[i - 1, ], w_bootparam_ef[i, ], r_oos_full, p),
                     calculate_to(w_boot_ef[i - 1, ], w_boot_ef[i, ], r_oos_full, p),
                     calculate_to(w_estim_pca_ef[i - 1, ], w_estim_pca_ef[i, ], r_oos_full, p),
                     calculate_to(w_factor_bootparam_pca_ef[i - 1, ], w_factor_bootparam_pca_ef[i, ], r_oos_full, p),
                     calculate_to(w_factor_boot_pca_ef[i - 1, ], w_factor_boot_pca_ef[i, ], r_oos_full, p),
                     calculate_to(w_estim_obs_ef[i - 1, ], w_estim_obs_ef[i, ], r_oos_full, p),
                     calculate_to(w_factor_bootparam_obs_ef[i - 1, ], w_factor_bootparam_obs_ef[i, ], r_oos_full, p),
                     calculate_to(w_factor_boot_obs_ef[i - 1, ], w_factor_boot_obs_ef[i, ], r_oos_full, p))
  }   
}

colnames(Rport) <- c("ew",
                     "mv_Markowitz", "mv_MichaudParam", "mv_MichaudNonP", 
                     "mv_MarkowitzPCA","mv_FactorParamPCA", "mv_FactorNonPPCA", 
                     "mv_MarkowitzObs","mv_FactorParamObs", "mv_FactorNonPObs",
                     "ew",
                     "ef_Markowitz", "ef_MichaudParam", "ef_MichaudNonP", 
                     "ef_MarkowitzPCA", "ef_FactorParamPCA", "ef_FactorNonPPCA", 
                     "ef_MarkowitzObs", "ef_FactorParamObs", "ef_FactorNonPObc")
write.table(Rport, paste0("Results/Rport_2_", InS, "_ibov.csv"), sep = ",")

colnames(to) <- colnames(Rport) 
write.table(to, paste0("Results/to_2_", InS, "_ibov.csv"), sep = ",")
 
colnames(sspw) <- colnames(Rport) 
write.table(sspw, paste0("Results/sspw_2_", InS, "_ibov.csv"), sep = ",")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#    Tables in Latex style     #  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Tables for short-selling MVP
Caption <- "Out-of-sample performance measures of the minimum variance portfolio with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[, 1:nmethods], 2, medidas, rf),
                     apply(to[, 1:nmethods], 2, mean),
                     apply(sspw[, 1:nmethods], 2, mean))
 row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
 colnames(oos_results) <- c("EW", "Markowitz", "Michaud Parametric", "Michaud Non-Parametric", "Markowitz-PCA", "Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Markowitz-Ibov", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
 t(oos_results) %>%
   knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
                table.envir = "table", label = "empirical_mvp") %>%
   save_kable(keep_tex = T, file = paste0("Results/mvp_", InS, "_ibov_5.tex"))

# Tables for short-selling EF
Caption <- "Out-of-sample performance measures of the eficient frontier with short-selling constraints: AV, SD, SR, ASR, SO, TO and SSPW stand for the average, standard deviation, Sharpe ratio, Adjusted Sharpe ratio, Sortino ratio, average turnover and average sum of squared portfolio weights, respectively."
oos_results <- rbind(apply(Rport[,(nmethods + 1):(2*nmethods)], 2, medidas, rf),  
                     apply(to[,(nmethods + 1):(2*nmethods)], 2, mean), 
                     apply(sspw[,(nmethods + 1):(2*nmethods)], 2, mean))
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO", "TO", "SSPW")
colnames(oos_results) <- c("EW","Markowitz", "Michaud Parametric", "Michaud Non-Parametric", "Markowitz-PCA", "Factor-Based Parametric PCA", "Factor-Based Non-Parametric PCA", "Markowitz-Ibov", "Factor-Based Parametric Ibov", "Factor-Based Non-Parametric Ibov")
t(oos_results) %>% 
  knitr::kable(digits = 4, format = "latex", align = "lccccccc", caption = Caption,
               table.envir = "table", label = "empirical_ef_2") %>% 
  save_kable(keep_tex = T, file = paste0("Results/ef_5_", InS, "_ibov.tex"))