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
pvalues_functions <- function(R, b) {
  mv_names <- colnames(R)[str_detect(colnames(R), "mv_")]
  tp_names <- colnames(R)[str_detect(colnames(R), "tp_")]
  ef_names <- colnames(R)[str_detect(colnames(R), "ef_")]
  pvalues_table <- matrix(NA, ncol = 6, nrow = (length(mv_names) - 1))
  colnames(pvalues_table) <- c("mv_sd", "mv_sr", "tp_sd", "tp_sr", "ef_sd", "ef_sr")
  row.names(pvalues_table) <- str_replace(mv_names[-1], "mv_", "")
  for (i in 1:(length(mv_names) - 1)) {
    Rtwo <- R %>% select("mv_Markowitz", mv_names[i + 1])
    pvalues_table[i,1] <- boot.time.inference.log.var(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
    pvalues_table[i,2] <- boot.time.inference(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
    Rtwo <- R %>% select("tp_Markowitz", tp_names[i + 1])
    pvalues_table[i,3] <- boot.time.inference.log.var(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
    pvalues_table[i,4] <- boot.time.inference(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
    Rtwo <- R %>% select("ef_Markowitz", ef_names[i + 1])
    pvalues_table[i,5] <- boot.time.inference.log.var(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
    pvalues_table[i,6] <- boot.time.inference(ret = Rtwo, b = 6, M = 5000, digits = 6)$p.Value
  }
  return(pvalues_table)
}

pvalues_functions_hac <- function(R, b) {
  mv_names <- colnames(R)[str_detect(colnames(R), "mv_")]
  tp_names <- colnames(R)[str_detect(colnames(R), "tp_")]
  ef_names <- colnames(R)[str_detect(colnames(R), "ef_")]
  pvalues_table <- matrix(NA, ncol = 6, nrow = (length(mv_names) - 1))
  colnames(pvalues_table) <- c("mv_sd", "mv_sr", "tp_sd", "tp_sr", "ef_sd", "ef_sr")
  row.names(pvalues_table) <- str_replace(mv_names[-1], "mv_", "")
  for (i in 1:(length(mv_names) - 1)) {
    Rtwo <- R %>% select("mv_Markowitz", mv_names[i + 1])
    pvalues_table[i,1] <- hac.inference.log.var(ret = Rtwo, digits = 6)$p.Values[2]
    pvalues_table[i,2] <- hac.inference(ret = Rtwo, digits = 6)$p.Values[2]
    Rtwo <- R %>% select("tp_Markowitz", tp_names[i + 1])
    pvalues_table[i,3] <- hac.inference.log.var(ret = Rtwo, digits = 6)$p.Values[2]
    pvalues_table[i,4] <- hac.inference(ret = Rtwo, digits = 6)$p.Values[2]
    Rtwo <- R %>% select("ef_Markowitz", ef_names[i + 1])
    pvalues_table[i,5] <- hac.inference.log.var(ret = Rtwo, digits = 6)$p.Values[2]
    pvalues_table[i,6] <- hac.inference(ret = Rtwo, digits = 6)$p.Values[2]
  }
  return(pvalues_table)
}

# ******************************************************

R60 = read.csv("./Results/Rport_60_ibov.csv")
R120 = read.csv("./Results/Rport_120_ibov.csv")
set.seed(321)
pvalues_ibov_60_hac <- pvalues_functions_hac(R60, b = 6)
pvalues_ibov_120_hac <- pvalues_functions_hac(R120, b = 6)

pvalues_ibov_60 <- pvalues_functions(R60, b = 6)
pvalues_ibov_120 <- pvalues_functions(R120, b = 6)

