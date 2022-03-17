#################################
###  Data Generating Process  ###
#################################

## Trying to mimic real monthly returns
setting_parameters <- function(p) {
  # Load Monthly returns
  monthly_data <- read_xlsx("Data/dados_mensais_IBRX.xlsx", na = "-")
  colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\r\ndo fechamento\r\nem 1 mês\r\nEm moeda orig\r\najust p/ prov\r\n", "")
  
  # Data Wrangling 
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
    filter(Data > '2012-12-31')
  
  monthly_data <- monthly_data %>% select(names(which(apply(is.na(monthly_data), 2, sum) == 0))) 
  monthly_data <- monthly_data[,2:(p + 1)]
  nobs <- nrow(monthly_data)
  monthly_data_c <- scale(monthly_data, center = TRUE, scale = FALSE)[,1:p]
  eigen_decomposition <- eigen(1 / (nobs - 1) * t(monthly_data_c) %*% monthly_data_c)
  n_factors <- POET::POETKhat(t(monthly_data_c))$K1HL
  #sum(eigen_decomposition$values > mean(eigen_decomposition$values))
  eigen_vectors <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  pca_factors <- as.matrix(monthly_data) %*% t(eigen_vectors)
  model <- lm(as.matrix(monthly_data[-1, ]) ~ pca_factors[-nobs, ])
  alpha <- coef(model)[1,]
  betas <- coef(model)[-1,]
  epsilon <- model$residuals
  Sigma_e <- 1/model$df.residual * t(epsilon) %*% epsilon
  while (any(eigen(Sigma_e)$values <= 1e-08)) {
    Sigma_e <- as.matrix(nearPD(Sigma_e, doSym = TRUE)$mat)
  }
  return(list(pca_factors, alpha, betas, Sigma_e))
}

# DGP
dgp <- function(nobs = 60, p = 10, tau = 100, set_par) {
  Sigma_e_pre <- set_par[[4]]
  betas_pre <- set_par[[3]]
  alpha_pre <- set_par[[2]]
  factors_pre <- set_par[[1]]
  n_oos <- 1000
  # Simulate epsilon
  Sigma_e <- rWishart(1, p, 1/p * Sigma_e_pre)[,,1]
  epsilon <- rmvnorm(nobs + n_oos , rep(0,p), Sigma_e)
  # Simulate alphas e betas
  alpha <- rmvnorm(1, alpha_pre, 1/tau * diag(p))
  betas <- rmvnorm(1, betas_pre, 1/tau * diag(p))
  # Simulate factor
  Sigma_factors <- var(factors_pre)*rchisq(1,1)
  k_factors <- matrix(rnorm(nobs + n_oos , mean = 0, sd = sqrt(Sigma_factors[1,1])), ncol = 1)
  # Simulated returns
  returns <- matrix(rep(alpha,nobs + n_oos), ncol = p, byrow = TRUE) + k_factors %*%  betas + epsilon
  mu <- alpha + mean(k_factors[c(1:nobs), ]) %*% beta
  # Simulated Covariance matrix
  Sigma <-  t(betas) %*% Sigma_factors %*% betas  + Sigma_e
  Sigma <- 0.5*Sigma + 0.5*t(Sigma)
  return(list(returns[1:nobs, ], Sigma, mu, returns[-c(1:nobs), ]))
}
