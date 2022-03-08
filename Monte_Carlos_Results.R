###################################
###     Monte Carlo Results     ###
###################################
library(stringr)
library(tidyr)
library(ggplot2)
setwd("./MC_d5_n60")

pattern <- "short_sales"
filenames <- list.files(pattern = "*.csv", full.names = TRUE)
filenames <- filenames[str_detect(filenames, pattern)]

w <- read.csv(filenames[1])[,-1]
n_methods <- length(filenames)
rmse <- matrix(NA, ncol = n_methods - 1, nrow = 1000)
mae <- matrix(NA, ncol = n_methods - 1, nrow = 1000)

for (j in 2:n_methods) {
  w_hat <- read.csv(filenames[j])[,-1]
  difer <- w - w_hat
  rmse[, j - 1] <- sqrt(apply(difer^2, 1, mean))
  mae[, j - 1] <- apply(abs(difer), 1, mean)
}

pattern2 <- paste0("./", pattern, "_mvp_w_")
names <- str_replace(str_replace(filenames[-1], paste0("5_60", ".csv"), ""), pattern2, "")
colnames(rmse) <- names
colnames(mae) <- names

apply(mae, 2, mean)
apply(rmse, 2, mean)

pivot_longer(data.frame(rmse), cols = 1:(n_methods - 1), names_to = "method") %>% 
  ggplot() + 
  geom_boxplot(aes(x = method, y = value, fill = method)) + 
  ylab("RMSE") + xlab("") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste0("rmse", pattern, ".pdf"), width = 30, height = 21, units = "cm")


pivot_longer(data.frame(mae), cols = 1:(n_methods - 1), names_to = "method") %>% 
  ggplot() + 
  geom_boxplot(aes(x = method, y = value, fill = method)) + 
  ylab("MAE") + xlab("") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(paste0("mae", pattern, ".pdf"), width = 30, height = 21, units = "cm")
