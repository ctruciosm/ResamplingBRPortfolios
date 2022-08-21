#################################
###   Out-of-sample figures   ###
#################################
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

to_60 <- read.csv("./Results/to_60_ibov.csv")
to_120 <- read.csv("./Results/to_60_ibov.csv")
sspw_60 <- read.csv("./Results/sspw_60_ibov.csv")
sspw_120 <- read.csv("./Results/sspw_120_ibov.csv")

sspw_60 <- sspw_60 %>% pivot_longer(cols = everything(), values_to = "sspw", names_to = "method") %>% mutate(size = "T = 60")
sspw_120 <- sspw_120 %>% pivot_longer(cols = everything(), values_to = "sspw", names_to = "method") %>% mutate(size = "T = 120")
sspw <- rbind(sspw_60, sspw_120) %>% 
  mutate(portfolio = str_sub(method, 1, 2), method = str_sub(method, 4)) %>% 
  mutate(portfolio = recode(portfolio,
                            "mv" = "Minimum Variance",
                            "tp" = "Tangency",
                            "ef" = "Mean-Variance")) %>% 
  mutate(method = recode(method,
                         "MichaudParam" = "Michaud P.",
                         "MichaudNonP" = "Michaud NP.",
                         "FactorParamPCA" = "Factor PCA P.",
                         "FactorNonPPCA" = "Factor PCA NP.",
                         "FactorParamObs" = "Factor Ibov P.",
                         "FactorNonPObs" = "Factor Ibov NP.",
                         "MarkowitzPCA" = "Markowitz PCA",
                         "MarkowitzObs" = "Markowitz Ibov")) %>% 
  mutate(method = factor(method, c("Markowitz", "Michaud P.", "Michaud NP.", "Markowitz PCA", "Factor PCA P.", "Factor PCA NP.", "Markowitz Ibov","Factor Ibov P.", "Factor Ibov NP.")))


to_60 <- to_60 %>% pivot_longer(cols = everything(), values_to = "to", names_to = "method") %>% mutate(size = "T = 60")
to_120 <- to_120 %>% pivot_longer(cols = everything(), values_to = "to", names_to = "method") %>% mutate(size = "T = 120")
to <- rbind(to_60, to_120) %>% 
  mutate(portfolio = str_sub(method, 1, 2), method = str_sub(method, 4)) %>% 
  mutate(portfolio = recode(portfolio,
                            "mv" = "Minimum Variance",
                            "tp" = "Tangency",
                            "ef" = "Mean-Variance")) %>% 
  mutate(method = recode(method,
                         "MichaudParam" = "Michaud P.",
                         "MichaudNonP" = "Michaud NP.",
                         "FactorParamPCA" = "Factor PCA P.",
                         "FactorNonPPCA" = "Factor PCA NP.",
                         "FactorParamObs" = "Factor Ibov P.",
                         "FactorNonPObs" = "Factor Ibov NP.",
                         "MarkowitzPCA" = "Markowitz PCA",
                         "MarkowitzObs" = "Markowitz Ibov")) %>% 
  mutate(method = factor(method, c("Markowitz", "Michaud P.", "Michaud NP.", "Markowitz PCA", "Factor PCA P.", "Factor PCA NP.", "Markowitz Ibov","Factor Ibov P.", "Factor Ibov NP.")))


ggplot(sspw) + 
  geom_boxplot(aes(x = factor(method), y = sspw, fill = factor(method))) + 
  facet_grid(size ~ portfolio) +
  ylab("SSPW") + xlab("") + 
  theme_bw() + 
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90))
ggsave("SSPW.pdf", height = 8.268, width = 11.693, limitsize = FALSE)



ggplot(to) + 
  geom_boxplot(aes(x = factor(method), y = to, fill = factor(method))) + 
  facet_grid(size ~ portfolio) +
  ylab("Turnover") + xlab("") + 
  theme_bw() + 
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90))
ggsave("to.pdf", height = 8.268, width = 11.693, limitsize = FALSE)


################################################
type <- "_tp"
w_boot_sd_60 <- read.csv("./Results/w_boot_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Michaud NP")
w_bootparam_sd_60 <- read.csv("./Results/w_bootparam_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Michaud P")
w_factor_boot_obs_sd_60 <- read.csv("./Results/w_factor_boot_obs_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Factor Obs NP")
w_factor_bootparam_obs_sd_60 <- read.csv("./Results/w_factor_bootparam_obs_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Factor Obs P")
w_factor_boot_pca_sd_60 <- read.csv("./Results/w_factor_boot_pca_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Factor PCA NP")
w_factor_bootparam_pca_sd_60 <- read.csv("./Results/w_factor_bootparam_pca_sd_120.csv") %>% select(ends_with(type)) %>% pivot_longer(cols = everything()) %>% mutate(name = str_replace(name, type,"")) %>% mutate(portfolio = "Factor PCA P")

weights <- data.frame(rbind(w_boot_sd_60, w_bootparam_sd_60, w_factor_boot_obs_sd_60, w_factor_bootparam_obs_sd_60, w_factor_boot_pca_sd_60, w_factor_bootparam_pca_sd_60))

ggplot(weights) + geom_boxplot(aes(x = portfolio, y = value, fill = portfolio)) +
  facet_wrap(.~name) + 
  xlab("") + ylab("") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90))
ggsave(paste0("weights_120", type, ".pdf"), height = 8.268, width = 11.693, limitsize = FALSE)
