#################################
###   Descriptive Statistics  ###
#################################
rm(list = ls())
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(lubridate)
library(kableExtra)
library(corrplot)
library(tseries)
library(ggplot2)
library(tsoutliers)

# Importing and Wrangling  Data
{
  monthly_data <- read_xlsx("Data/economatica_ibov.xlsx", na = "-")
  colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\ndo fechamento\nem 1 mês\nEm moeda orig\najust p/ prov\n", "")
  
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
    filter(Data >= '2000-01-01', Data <= '2022-02-01')
  
  monthly_data_dates <- monthly_data %>% 
    select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) 
  
  ibovespa <- monthly_data_dates %>% select(IBOV)
  monthly_data <- monthly_data_dates %>% select(-IBOV)
}

monthly_data_longer <- pivot_longer(monthly_data, cols = ALPA4:VALE3, values_to = "returns", names_to = "assets")

assets_names <- 
  monthly_data_longer %>% 
  group_by(assets) %>% 
  summarise(LjungBox = Box.test(returns, type =  "Ljung-Box")$p.value,
            LjungBox2 = Box.test(returns^2, type =  "Ljung-Box")$p.value) %>% 
  filter(LjungBox > 0.05, LjungBox2 > 0.05) %>% 
  select(assets)

monthly_data <-  monthly_data_dates %>% select(Data, as.vector(as.matrix(assets_names)))
monthly_data_longer <- pivot_longer(monthly_data, cols = ALPA4:VIVT3, values_to = "returns", names_to = "assets")


# Figure 1
monthly_data_longer %>% 
  ggplot() + geom_line(aes(x = Data, y = returns, color = assets)) + 
  xlab(" ") + 
  ylab("Monthly returns") + theme_bw() + theme(legend.position = "none", legend.title = element_blank())
ggsave("ibov_monthly_returns.pdf", width = 37, height = 19, units = "cm")           


# Table 1
monthly_data_longer %>% 
  select(-Data) %>% 
  group_by(assets) %>% 
  summarise(minimum = min(returns),
            maximum = max(returns),
            mean = mean(returns),
            sd = sd(returns),
            skewness = moments::skewness(returns),
            kurtosis = moments::kurtosis(returns),
            JB = jarque.bera.test(returns)$p.value) %>% 
  knitr::kable(digits = 3, format = "latex", align = "lccccccc",
               table.envir = "table", label = "descriptive_statistics") %>% 
  save_kable(keep_tex = T, file = paste0("Results/descriptive_statistics_ibov.tex"))