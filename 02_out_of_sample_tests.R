#################################
###    Out-of-sample tests    ###
#################################
library(dplyr)
load("Results/Var/Var.RData")
load("Results/SharpeR/Sharpe.RData")

block <- 12

Ropt_ls_0 <- read.csv("Results/Rport_2_120_ls_0_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_1 <- read.csv("Results/Rport_2_120_ls_1_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_2 <- read.csv("Results/Rport_2_120_ls_2_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_3 <- read.csv("Results/Rport_2_120_ls_3_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_4 <- read.csv("Results/Rport_2_120_ls_4_ibrx.csv")[, c(2:10, 1)]

# Minimum Variance
p_values_mv <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_mv[i - 1, 1] <- boot.time.inference.log.var(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_mv[i - 1, 2] <- boot.time.inference.log.var(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_mv[i - 1, 3] <- boot.time.inference.log.var(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_mv[i - 1, 4] <- boot.time.inference.log.var(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_mv[i - 1, 5] <- boot.time.inference.log.var(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_mv) <- colnames(Ropt_ls_0)[-1]

Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_mv_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_mv_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}

# Lambda 2
Ropt_ls_0 <- read.csv("Results/Rport_2_120_ls_0_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_1 <- read.csv("Results/Rport_2_120_ls_1_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_2 <- read.csv("Results/Rport_2_120_ls_2_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_3 <- read.csv("Results/Rport_2_120_ls_3_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_4 <- read.csv("Results/Rport_2_120_ls_4_ibrx.csv")[, c(12:20, 11)]

p_values_2 <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_2[i - 1, 1] <- boot.time.inference(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_2[i - 1, 2] <- boot.time.inference(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_2[i - 1, 3] <- boot.time.inference(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_2[i - 1, 4] <- boot.time.inference(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_2[i - 1, 5] <- boot.time.inference(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_2) <- colnames(Ropt_ls_0)[-1]

Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_2_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_2_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}


# lambda 5
Ropt_ls_0 <- read.csv("Results/Rport_5_7_120_ls_0_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_1 <- read.csv("Results/Rport_5_7_120_ls_1_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_2 <- read.csv("Results/Rport_5_7_120_ls_2_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_3 <- read.csv("Results/Rport_5_7_120_ls_3_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_4 <- read.csv("Results/Rport_5_7_120_ls_4_ibrx.csv")[, c(2:10, 1)]

p_values_5 <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_5[i - 1, 1] <- boot.time.inference(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_5[i - 1, 2] <- boot.time.inference(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_5[i - 1, 3] <- boot.time.inference(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_5[i - 1, 4] <- boot.time.inference(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_5[i - 1, 5] <- boot.time.inference(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_5) <- colnames(Ropt_ls_0)[-1]

Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_5_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_5_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}



# lambda 7
Ropt_ls_0 <- read.csv("Results/Rport_5_7_120_ls_0_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_1 <- read.csv("Results/Rport_5_7_120_ls_1_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_2 <- read.csv("Results/Rport_5_7_120_ls_2_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_3 <- read.csv("Results/Rport_5_7_120_ls_3_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_4 <- read.csv("Results/Rport_5_7_120_ls_4_ibrx.csv")[, c(12:20, 11)]

# Tests within panel (panel = cov estimator)
p_values_7 <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_7[i - 1, 1] <- boot.time.inference(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_7[i - 1, 2] <- boot.time.inference(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_7[i - 1, 3] <- boot.time.inference(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_7[i - 1, 4] <- boot.time.inference(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_7[i - 1, 5] <- boot.time.inference(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_7) <- colnames(Ropt_ls_0)[-1]


Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_7_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_7_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}


# lambda 05
Ropt_ls_0 <- read.csv("Results/Rport_05_1_120_ls_0_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_1 <- read.csv("Results/Rport_05_1_120_ls_1_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_2 <- read.csv("Results/Rport_05_1_120_ls_2_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_3 <- read.csv("Results/Rport_05_1_120_ls_3_ibrx.csv")[, c(2:10, 1)]
Ropt_ls_4 <- read.csv("Results/Rport_05_1_120_ls_4_ibrx.csv")[, c(2:10, 1)]

# Tests within panel (panel = cov estimator)
p_values_05 <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_05[i - 1, 1] <- boot.time.inference(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_05[i - 1, 2] <- boot.time.inference(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_05[i - 1, 3] <- boot.time.inference(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_05[i - 1, 4] <- boot.time.inference(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_05[i - 1, 5] <- boot.time.inference(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_05) <- colnames(Ropt_ls_0)[-1]


Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_05_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_05_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}


# lambda 11
Ropt_ls_0 <- read.csv("Results/Rport_05_1_120_ls_0_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_1 <- read.csv("Results/Rport_05_1_120_ls_1_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_2 <- read.csv("Results/Rport_05_1_120_ls_2_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_3 <- read.csv("Results/Rport_05_1_120_ls_3_ibrx.csv")[, c(12:20, 11)]
Ropt_ls_4 <- read.csv("Results/Rport_05_1_120_ls_4_ibrx.csv")[, c(12:20, 11)]

# Tests within panel (panel = cov estimator)
p_values_1 <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:10) {
  set.seed(i)
  p_values_1[i - 1, 1] <- boot.time.inference(Ropt_ls_0[, c(1, i)], b = block, M = 500)$p.Value
  p_values_1[i - 1, 2] <- boot.time.inference(Ropt_ls_1[, c(1, i)], b = block, M = 500)$p.Value
  p_values_1[i - 1, 3] <- boot.time.inference(Ropt_ls_2[, c(1, i)], b = block, M = 500)$p.Value
  p_values_1[i - 1, 4] <- boot.time.inference(Ropt_ls_3[, c(1, i)], b = block, M = 500)$p.Value
  p_values_1[i - 1, 5] <- boot.time.inference(Ropt_ls_4[, c(1, i)], b = block, M = 500)$p.Value
}
rownames(p_values_1) <- colnames(Ropt_ls_0)[-1]


Ropt_mv <- cbind(Ropt_ls_0[, 1:9], Ropt_ls_1[, 1:9], Ropt_ls_2[, 1:9], Ropt_ls_3[, 1:9], Ropt_ls_4)
p_values_1_full <- c()
for (i in 2:46) {
  set.seed(i)
  p_values_1_full[i - 1] <- boot.time.inference.log.var(Ropt_mv[, c(1, i)], b = block, M = 500)$p.Value
}





p_values_mv
p_values_2
p_values_5
p_values_7
p_values_05
p_values_1

p_values_mv_full
p_values_2_full
p_values_5_full
p_values_7_full
p_values_05_full
p_values_1_full
