#################################
###    Auxiliary Functions    ###
#################################

covariance_method <- function(x, method = 1) {
  Sigma_est <- switch(method,
                      cov(x),
                      CovMcd(x)$cov,
                      CovMest(x)$cov,
                      CovMrcd(x)$cov,
                      CovMve(x)$cov,
                      CovOgk(x)$cov,
                      CovSest(x)$cov,
                      nlShrinkLWEst(x),
                      nlshrink_cov(x),
                      covEstimation(x, control = list(type = 'lw')),
                      covEstimation(x, control = list(type = 'oneparm')),
                      covEstimation(x, control = list(type = 'const')),
                      covEstimation(x, control = list(type = 'cor')),
                      covEstimation(x, control = list(type = 'diag')),
                      covEstimation(x, control = list(type = 'large')))
  return(Sigma_est)
}


calculate_portfolio_weights <- function(x, constrains_opt, nboot) {
  w_estim <- optimalPortfolio(Sigma = cov(x), mu = colMeans(x), control = constrains_opt)
  w_bootparam <- michaud_parametric_bootstrap(x, B = nboot, option_list = constrains_opt) 
  w_boot <- michaud_bootstrap(x, B = nboot, option_list = constrains_opt) 
  k <- POETKhat(t(x))$K1HL
  w_factor_bootparam <- factor_parametric_bootstrap(x, B = nboot, n_factors  = k, option_list = constrains_opt) 
  w_factor_boot <- factor_bootstrap(x, B = nboot, n_factors  = k, option_list = constrains_opt) 
  return(rbind(w_estim, w_bootparam, w_boot, w_factor_bootparam, w_factor_boot))
}