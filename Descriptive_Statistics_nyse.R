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

# Importing and Wrangling  Data
{
  monthly_data <- read_xlsx("Data/economatica_nyse.xlsx",  skip = 3)
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
    mutate(Date = lubridate::my(Data)) %>% 
    mutate_if(is.character, as.numeric) %>% 
    filter(Date >= '1995-01-01', Date <= '2022-01-01')
  monthly_data <- monthly_data %>% 
    select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) %>% 
    filter(Date >= '1995-01-01', Date <= '2022-01-01') %>% 
    select(-Date)
}

monthly_data_longer <- pivot_longer(monthly_data, cols = everything(), values_to = "returns", names_to = "assets")

# Figure 1
monthly_data_longer %>% 
  mutate(Dates = rep(monthly_data_dates$Data[-1], ncol(monthly_data))) %>%
  ggplot() + geom_line(aes(x = Dates, y = returns, color = assets)) + 
  ylab("Monthly returns") + theme(legend.position = "bottom", legend.title = element_blank())
           
monthly_data_longer %>% 
  mutate(Dates = rep(monthly_data_dates$Data[-1], ncol(monthly_data))) %>%
  ggplot() + geom_line(aes(x = Dates, y = returns, color = assets)) + 
  facet_wrap(.~ assets) + 
  ylab("Monthly returns") + theme(legend.position = "none")

# Table 1
monthly_data_longer %>% 
  group_by(assets) %>% 
  summarise(minimum = min(returns),
            maximum = max(returns),
            mean = mean(returns),
            sd = sd(returns),
            skewness = moments::skewness(returns),
            kurtosis = moments::kurtosis(returns),
            Shapiro = shapiro.test(returns)$p.value,
            LjungBox = Box.test(returns, lag = 1, type =  "Ljung-Box")$p.value,
            LjungBox2 = Box.test(returns^2, lag = 1, type =  "Ljung-Box")$p.value) %>% 
  knitr::kable(digits = 3, format = "latex", align = "lccccccccc",
               table.envir = "table", label = "descriptive_statistics") %>% 
  save_kable(keep_tex = T, file = paste0("Results/descriptive_statistics.tex"))

# Figure Appendix
corrplot(cor(monthly_data), type = 'lower', number.cex = 0.6, addCoef.col = 'grey50', diag = FALSE, method = 'color')
