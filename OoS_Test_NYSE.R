##################################
###      Bootstrap Tests       ###
##################################
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
load("./Results/Var/Var.RData")
load("./Results/SharpeR/Sharpe.RData")
# ******************************************************
pvalues_functions <- function(R) {
  ssc_names <- colnames(R)[str_detect(colnames(R), "ssc_")]
  luc_names <- colnames(R)[str_detect(colnames(R), "luc_")]
  pvalues_table <- matrix(NA, ncol = 4, nrow = (length(luc_names) - 1))
  colnames(pvalues_table) <- c("ssc_sd", "ssc_sr", "luc_sd", "luc_sr")
  row.names(pvalues_table) <- str_replace(ssc_names[-1], "ssc_", "")
  for (i in 1:(length(luc_names) - 1)) {
    Rtwo <- R %>% select("ssc_Markowitz", ssc_names[i + 1])
    pvalues_table[i,1] <- boot.time.inference.log.var(ret = Rtwo, b = 10, M = 5000)$p.Value
    pvalues_table[i,2] <- boot.time.inference(ret = Rtwo, b = 10, M = 5000)$p.Value
    Rtwo <- R %>% select("luc_Markowitz", luc_names[i + 1])
    pvalues_table[i,3] <- boot.time.inference.log.var(ret = Rtwo, b = 10, M = 5000)$p.Value
    pvalues_table[i,4] <- boot.time.inference(ret = Rtwo, b = 10, M = 5000)$p.Value
  }
  return(pvalues_table)
}
# ******************************************************

R60 = read.csv("./Results/Rport_60_nyse.csv")
R120 = read.csv("./Results/Rport_120_nyse.csv")
set.seed(123)
pvalues_nyse_60 <- pvalues_functions(R60)
pvalues_nyse_120 <- pvalues_functions(R120)



R60 = read.csv("./Results/Rport_60_ibov.csv")
R120 = read.csv("./Results/Rport_120_ibov.csv")
set.seed(246)
pvalues_ibov_60 <- pvalues_functions(R60)
pvalues_ibov_120 <- pvalues_functions(R120)






R %>% 
  select(starts_with("ssc_")) %>% 
  mutate_if(is.numeric, cumsum) %>% 
  mutate(Index = 1:260) %>%
  select(Index, everything()) %>% 
  pivot_longer(cols = ssc_Markowitz:ssc_FactorNonPFF5f, names_to = "Methods", values_to = "cumreturns") %>% 
  ggplot() +
  geom_line(aes(x = Index, y = cumreturns, color = Methods))
