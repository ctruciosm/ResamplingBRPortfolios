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
                         "FactorParamFF" = "Factor Ibov P.",
                         "FactorNonPFF" = "Factor Ibov NP.")) %>% 
  mutate(method = factor(method, c("Markowitz", "Michaud P.", "Michaud NP.", "Factor PCA P.", "Factor PCA NP.", "Factor Ibov P.", "Factor Ibov NP.")))


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
                         "FactorParamFF" = "Factor Ibov P.",
                         "FactorNonPFF" = "Factor Ibov NP.")) %>% 
  mutate(method = factor(method, c("Markowitz", "Michaud P.", "Michaud NP.", "Factor PCA P.", "Factor PCA NP.", "Factor Ibov P.", "Factor Ibov NP.")))



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