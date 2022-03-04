###################################
###    Auxiliary Functions      ###
###################################

# Setting Parameters
library(tidyr)
library(readxl)
library(stringr)
library(dplyr)
## Similar values to the observed empirically
setting_parameters <- function(p) {
  # Load Monthly returns
  monthly_data <- read_xlsx("Data/dados_mensais_IBRX.xlsx", na = "-")
  colnames(monthly_data) <- str_replace(colnames(monthly_data), "Retorno\r\ndo fechamento\r\nem 1 mês\r\nEm moeda orig\r\najust p/ prov\r\n", "")
  
  # Data Wrangling 
  monthly_data <- monthly_data %>%
    mutate(Data = str_replace(Data, "Jan", "01"), Data = str_replace(Data, "Fev", "02"),
           Data = str_replace(Data, "Mar", "03"), Data = str_replace(Data, "Abr", "04"),
           Data = str_replace(Data, "Mai", "05"), Data = str_replace(Data, "Jun", "06"),
           Data = str_replace(Data, "Jul", "07"), Data = str_replace(Data, "Ago", "08"),
           Data = str_replace(Data, "Set", "09"), Data = str_replace(Data, "Out", "10"),
           Data = str_replace(Data, "Nov", "11"), Data = str_replace(Data, "Dez", "12")) %>% 
    mutate(Data = lubridate::my(Data)) %>% 
    filter(Data > '2012-12-31')
  
  monthly_data <- monthly_data %>% select(names(which(apply(is.na(monthly_data), 2, sum) == 0)))
  monthly_data <- monthly_data[,sample(2:ncol(monthly_data), p, replace = FALSE)]
  monthly_data_centred <- scale(monthly_data, center = TRUE, scale = FALSE)
  pca <- princomp(monthly_data_centred, cor = FALSE)
  # Bai and Ng as well as Hallin and Liska point out 1 factor
  onefactor <- matrix(pca$scores[,1], ncol = 1) 
  betas <- matrix(pca$loadings[,1],nrow = 1)
  model <- lm(as.matrix(monthly_data_centred) ~ onefactor-1)
  epsilon <- residuals(model)
  #epsilon <-  monthly_data_centred - onefactor %*% betas
  Sigma_e <- cov(epsilon)
  return(list(onefactor, betas, Sigma_e))
}

# DGP
dgp <- function(n = 60, p = 10, includemean = FALSE, tau = 100, set_par) {
  Sigma_e_pre <- sp[[3]]
  betas_pre <- sp[[2]]
  onefactor_pre <- sp[[1]]
  # Simulate epsilon
  Sigma_e <- rWishart(1, p, 1/p * Sigma_e_pre)[,,1]
  epsilon <- rmvnorm(n, rep(0,p), Sigma_e)
  # Simulate betas
  betas <- rmvnorm(1, betas_pre, 1/tau * diag(p))
  betas <- matrix(betas/sqrt(sum(betas^2)), nrow = 1)
  # Simulate factor
  sigma_factor <- sd(onefactor_pre)*rchisq(1,1)
  onefactor <- matrix(rnorm(n, mean = 0, sd = sigma_factor), ncol = 1)
  # Simulated returns
  returns <- onefactor %*%  betas + epsilon
  mu <- rep(0, p)
  if (includemean) {
    mu <- matrix(rnorm(p, mean = 0, sd = 1/p), ncol = p, nrow = n, byrow = TRUE)
    returns <- mu + returns
  } 
  # Simulated Covariance matrix
  Sigma <- sigma_factor * t(betas) %*% betas  + Sigma_e
  return(list(returns, Sigma, mu))
}

# Michaud Parametric Bootstrap
michaud_parametric_bootstrap <- function(x, B = 1000, p) {
  mu_hat <- apply(x, 2, mean)
  Sigma_hat <- cov(x)
  n <- nrow(x)

  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- rmvnorm(n, mu_hat, Sigma_hat)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Michaud Bootstrap
michaud_bootstrap <- function(x, B = 1000, p) {
  mu_hat <- apply(x, 2, mean)
  Sigma_hat <- cov(x)
  n <- nrow(x)
  eta_hat <- scale(x, center = TRUE, scale = FALSE) %*% solve(chol(Sigma_hat))
  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  
  for (j in 1:B) {
    eta_boot <- eta_hat[sample(1:n,n, replace = TRUE),]
    returns_boot <- mu_hat + eta_boot %*% chol(Sigma_hat)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat,  control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Conditional Factor Parametric Bootstrap
factor_parametric_bootstrap <- function(x, B = 1000, p) {
  n <- nrow(x)
  pca <- princomp(x)
  pca_factors <- pca$scores[,1]
  linear_model <- lm(as.matrix(x) ~ pca_factors)
  betas <- coef(linear_model)
  epsilon <- residuals(linear_model) 
  Sigma_e <- cov(epsilon)

  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    epsilon_boot <- rmvnorm(n, rep(0,p), Sigma_e)
    returns_boot <- cbind(rep(1,n),pca_factors) %*% betas + epsilon_boot
    #betas[1,] + mean(pca_factors)*betas[2,] + apply(epsilon_boot,2,mean)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    #var(pca_factors) * t(matrix(betas[2,], nrow = 1)) %*% matrix(betas[2,], nrow = 1)  + cov(epsilon_boot)
    #
    #t(matrix(betas, nrow = 2)) %*% cov(cbind(rep(1,n),pca_factors)) %*% matrix(betas, nrow = 2)  + cov(epsilon_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat, control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}

# Conditional Factor Bootstrap
factor_bootstrap <- function(x, B = 1000, p) {
  n <- nrow(x)
  pca <- princomp(x)
  pca_factors <- pca$scores[,1]
  linear_model <- lm(as.matrix(x) ~ pca_factors)
  betas <- coef(linear_model)
  epsilon <- residuals(linear_model) 
  Sigma_e <- cov(epsilon)
  eta <- epsilon %*% solve(chol(Sigma_e))

  umvw_boot <- matrix(NA, ncol = p, nrow = B)
  for (j in 1:B) {
    returns_boot <- cbind(rep(1,n),pca_factors) %*% betas + eta[sample(1:n, n, replace = TRUE),] %*% chol(Sigma_e)
    mu_boot_hat <- apply(returns_boot, 2, mean)
    Sigma_boot_hat <- cov(returns_boot)
    umvw_boot[j, ] <- optimalPortfolio(Sigma = Sigma_boot_hat, mu = mu_boot_hat,  control = list(type = 'minvol'))
  }
  return(w = apply(umvw_boot, 2, mean))
}


