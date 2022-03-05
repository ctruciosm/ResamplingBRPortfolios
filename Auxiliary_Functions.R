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
