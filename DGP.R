#################################
###  Data Generating Process  ###
#################################

## Trying to mimic real monthly returns
setting_parameters <- function(p) {
  # Load Monthly returns
  monthly_data <- read_xlsx("Data/economatica.xlsx", na = "-")
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
    filter(Data >= '2012-01-01', Data <= '2022-01-01')
    monthly_data <- monthly_data %>% select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) 
    ibovespa <- monthly_data %>% select(IBOV)
    monthly_data <- monthly_data %>% select(-Data, -IBOV)
    
    model <- lm(as.matrix(monthly_data[-1, 1:p]) ~ as.matrix(ibovespa[-nrow(monthly_data), ]))
    alpha_hat <- coef(model)[1,]
    beta_hat <- matrix(coef(model)[-1, ], ncol = p)
    Sigma_e_hat <- (1/model$df.residual) * t(model$residuals) %*% model$residuals
    mu_F = mean(unlist(ibovespa[-nrow(monthly_data),]))
    sigma2_F = var(unlist(ibovespa[-nrow(monthly_data),]))
  
  return(list(alpha_hat, beta_hat, Sigma_e_hat, mu_F, sigma2_F))
}

# DGP
dgp <- function(nobs = 60, p = 10, set_par) {
  alpha_hat <- set_par[[1]]
  beta_hat <- set_par[[2]]
  Sigma_e_hat <- set_par[[3]]
  mu_F_hat <- set_par[[4]]
  sigma2_F_hat <- set_par[[5]]
  n_oos <- 1000
  # Simulate epsilon
  Sigma_e <- rWishart(1, p, 1/p * Sigma_e_hat)[,,1]
  epsilon <- rmvnorm(nobs + n_oos , rep(0,p), Sigma_e)
  # Simulate alphas e betas
  alpha <- rmvnorm(1, alpha_hat, diag(p))
  betas <- rmvnorm(1, beta_hat, diag(p))
  # Simulate factor
  Sigma2_F <- sigma2_F_hat*rchisq(1,1)
  mu_f <- rnorm(1,mu_F_hat, 1)
  one_factor <- matrix(rnorm(nobs + n_oos , mean = mu_f, sd = sqrt(Sigma2_F)), ncol = 1)
  # Simulated returns
  returns <- matrix(rep(alpha,nobs + n_oos), ncol = p, byrow = TRUE) + one_factor %*%  betas + epsilon
  mu <- alpha + mu_f %*% betas
  Sigma <-  Sigma2_F*t(betas) %*% betas  + Sigma_e
  Sigma <- 0.5*Sigma + 0.5*t(Sigma)
  return(list(returns[1:nobs, ], Sigma, mu, returns[-c(1:nobs), ], c(one_factor[2:nobs, ],0)))
}
