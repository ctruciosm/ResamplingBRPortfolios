#################################
###   Empirical Application   ###
#################################

library(dplyr)
source("Resampling_Techniques.R")
source("Performance_Measures.R")

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
constrains_opt <- list(type = 'minvol', constraint = 'lo')
w_estim = w_bootparam = w_boot = w_factor_bootparam = w_factor_boot = matrix(NA, ncol = p, nrow = OoS)
Rport <- matrix(NA, ncol = 5, nrow = OoS)
for (i in 1:OoS) {
  print(i)
  set.seed(i + 123)
  r <- as.matrix(monthly_data[i:(InS - 1 + i), ])
  w_estim[i,] <- optimalPortfolio(Sigma = cov(r), mu = colMeans(r), control = constrains_opt)
  w_bootparam[i,] <- michaud_parametric_bootstrap(r, B = nboot, option_list = constrains_opt) 
  w_boot[i,] <- michaud_bootstrap(r, B = nboot, option_list = constrains_opt) 
  k <- POET::POETKhat(t(r))$K1HL
  w_factor_bootparam[i,] <- factor_parametric_bootstrap(r, B = nboot, n_factors  = k, option_list = constrains_opt) 
  w_factor_boot[i,] <- factor_bootstrap(r, B = nboot, n_factors  = k, option_list = constrains_opt) 
  OoS_returns <- as.numeric(monthly_data[InS - 1 + i + 1, ]) 
  Rport[i,] <- c(OoS_returns %*% w_estim[i,], OoS_returns %*% w_bootparam[i,], 
                OoS_returns %*% w_boot[i,], OoS_returns %*% w_factor_bootparam[i,],
                OoS_returns %*% w_factor_boot[i,])
  
}

oos_results <- apply(Rport,2,medidas)
row.names(oos_results) <- c("AV", "SD", "SR", "ASR", "SO")
colnames(oos_results) <- c("Markowitz", "Michaud Parametric", "Michaud Non-Parametric",
                            "Factor-Based Parametric", "Factor-Based Non-Parametric")

xtable::xtable(t(oos_results),digits = 4)




