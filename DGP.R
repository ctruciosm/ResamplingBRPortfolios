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
  monthly_data_centred <- scale(monthly_data, center = TRUE, scale = FALSE)[,1:p]
  eigen_decomposition <- eigen(cov(monthly_data_centred))
  n_factors <- POET::POETKhat(t(monthly_data_centred))$K1HL
  #sum(eigen_decomposition$values > mean(eigen_decomposition$values))
  betas <- matrix(t(eigen_decomposition$vectors[,1:n_factors]), nrow = n_factors)
  pca_factors <- monthly_data_centred %*% t(betas)
  epsilon <-  monthly_data_centred - pca_factors %*% betas
  Sigma_e <- cov(epsilon)
  while (any(eigen(Sigma_e)$values <= 1e-08)) {
    Sigma_e <- as.matrix(nearPD(Sigma_e, doSym = TRUE)$mat)
  }
  return(list(pca_factors, betas, Sigma_e))
}

# DGP
dgp <- function(nobs = 60, p = 10, includemean = FALSE, tau = 100, set_par) {
  Sigma_e_pre <- set_par[[3]]
  betas_pre <- set_par[[2]]
  factors_pre <- set_par[[1]]
  k <- ncol(factors_pre)
  # Simulate epsilon
  Sigma_e <- rWishart(1, p, 1/p * Sigma_e_pre)[,,1]
  epsilon <- rmvnorm(nobs, rep(0,p), Sigma_e)
  # Simulate betas
  betas <- rmvnorm(1, betas_pre, 1/tau * diag(p))
  betas <- matrix(betas/sqrt(sum(betas^2)), nrow = k)
  # Simulate factor
  if (k == 1) {
    Sigma_factors <- var(factors_pre)*rchisq(1,1)
    k_factors <- matrix(rnorm(nobs, mean = 0, sd = sqrt(Sigma_factors[1,1])), ncol = k)
  } else{
    Sigma_factors <- rWishart(1, k, 1/k * cov(factors_pre))[,,1]
    k_factors <- rmvnorm(nobs, rep(0,p), Sigma_factors)
  }
  # Simulated returns
  returns <- k_factors %*%  betas + epsilon
  mu <- rep(0, p)
  if (includemean) {
    mu <- matrix(rnorm(p, mean = 0, sd = 1/p), ncol = p, nrow = nobs, byrow = TRUE)
    returns <- mu + returns
  } 
  # Simulated Covariance matrix
  Sigma <-  t(betas) %*% Sigma_factors %*% betas  + Sigma_e
  Sigma <- 0.5*Sigma + 0.5*t(Sigma)
  return(list(returns, Sigma, mu))
}
