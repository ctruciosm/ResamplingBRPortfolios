################################
###     Data Wrangling       ###
################################

library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(lubridate)
library(nortsTest)


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
  dplyr::filter(Data >= '2000-01-01', Data <= '2022-02-01')


monthly_data_dates <- monthly_data %>% 
  dplyr::select(names(which(apply(is.na(monthly_data), 2, sum) == 0)))

monthly_data_dates %>% dplyr::select(IBOV) %>% write.table("ibovespa.csv", sep = ",")
monthly_data <- monthly_data_dates %>% dplyr::select(-Data, -IBOV)

monthly_data_longer <- pivot_longer(monthly_data, cols = everything(), values_to = "returns", names_to = "assets")

assets_names <- 
  monthly_data_longer1 %>% 
  group_by(assets) %>% 
  summarise(LjungBox = Box.test(returns, type =  "Ljung-Box")$p.value,
            LjungBox2 = Box.test(returns^2, type =  "Ljung-Box")$p.value) %>% 
  dplyr::filter(LjungBox > 0.05, LjungBox2 > 0.05) %>% 
  dplyr::select(assets)

monthly_data %>% 
  dplyr::select(as.vector(as.matrix(assets_names))) %>% 
  write.table("monthly_data.csv", sep = ",")
