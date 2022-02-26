###################################
###      Simulation Setup       ###
###################################
library(mvtnorm)
library(RiskPortfolios)

n_assets <- 10
t_assets <- 100
mu <- rnorm(10)
H <- rWishart(1, 20, toeplitz((10:1)/10))[,,1]
nu <- 26
tau <- 5




H_sim <- rWishart(1, nu, 1/nu * H)[,,1]
mu_sim <- rmvnorm(1, mu, 1/tau * H_sim)
for (i in 1:100) {
  set.seed(i+123)
  returns_sim <- rmvnorm(t_assets, mu_sim, H_sim)
  mu_sim_hat <- apply(returns_sim, 2, mean)
  H_sim_hat <- cov(returns_sim)
  # Markowitz weights
  umvw <- optimalPortfolio(Sigma = H_sim_hat, mu = mu_sim_hat, control = list(type = 'mv'))
  # Michaud weights
  umvw_michaud_boot <- matrix(NA, ncol = n_assets, nrow = 500)
  for (j in 1:500) {
    returns_sim_boot <- rmvnorm(ntotal, mu_sim_hat, H_sim_hat)
    mu_sim_boot_hat <- apply(returns_sim_boot, 2, mean)
    H_sim_boot_hat <- cov(returns_sim_boot)
    umvw_michaud_boot[j, ] <- optimalPortfolio(Sigma = H_sim_boot_hat, mu = mu_sim_boot_hat, control = list(type = 'mv'))
  }
  weights_michaud_boot <- apply(umvw_michaud_boot, 2, mean)
}
