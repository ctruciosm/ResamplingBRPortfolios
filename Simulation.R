###################################
###      Simulation Setup       ###
###################################
library(mvtnorm)

mu <- rnorm(10)
H <- rWishart(1, 20, toeplitz((10:1)/10))[,,1]
nu <- 26
tau <- 5
ntotal <- 100



H_sim <- rWishart(1, nu, 1/nu * H)[,,1]
mu_sim <- rmvnorm(1, mu, 1/tau * H_sim)
for (i in 1:100) {
  returns_sim <- rmvnorm(ntotal, mu_sim, H_sim)
  
}
