#################################
###    Auxiliary Functions    ###
#################################


tangency.portfolio <-
  function(er,cov.mat,risk.free, shorts=TRUE)
  {
    call <- match.call()
    
    #
    # check for valid inputs
    #
    asset.names <- names(er)
    if(risk.free < 0)
      stop("Risk-free rate must be positive")
    er <- as.vector(er)
    cov.mat <- as.matrix(cov.mat)
    N <- length(er)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum variance portfolio
    #
    gmin.port <- globalMin.portfolio(er, cov.mat, shorts=shorts)
    if(gmin.port$er < risk.free)
      stop("Risk-free rate greater than avg return on global minimum variance portfolio")
    
    # 
    # compute tangency portfolio
    #
    if(shorts==TRUE){
      cov.mat.inv <- solve(cov.mat)
      w.t <- cov.mat.inv %*% (er - risk.free) # tangency portfolio
      w.t <- as.vector(w.t/sum(w.t))          # normalize weights
    } else if(shorts==FALSE){
      Dmat <- 2*cov.mat
      dvec <- rep.int(0, N)
      er.excess <- er - risk.free
      Amat <- cbind(er.excess, diag(1,N))
      bvec <- c(1, rep(0,N))
      result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w.t <- round(result$solution/sum(result$solution), 6)
    } else {
      stop("Shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    names(w.t) <- asset.names
    er.t <- crossprod(w.t,er)
    sd.t <- sqrt(t(w.t) %*% cov.mat %*% w.t)
    tan.port <- list("call" = call,
                     "er" = as.vector(er.t),
                     "sd" = as.vector(sd.t),
                     "weights" = w.t)
    class(tan.port) <- "portfolio"
    return(tan.port)
  }




globalMin.portfolio <-
  function(er, cov.mat, shorts=TRUE)
  {
    call <- match.call()
    
    #
    # check for valid inputs
    #
    asset.names <- names(er)
    er <- as.vector(er) # assign names if none exist
    cov.mat <- as.matrix(cov.mat)
    N <- length(er)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum portfolio
    #
    if(shorts==TRUE){
      cov.mat.inv <- solve(cov.mat)
      one.vec <- rep(1,N)
      w.gmin <- rowSums(cov.mat.inv) / sum(cov.mat.inv)
      w.gmin <- as.vector(w.gmin)
    } else if(shorts==FALSE){
      Dmat <- 2*cov.mat
      dvec <- rep.int(0, N)
      Amat <- cbind(rep(1,N), diag(1,N))
      bvec <- c(1, rep(0,N))
      result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w.gmin <- round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    names(w.gmin) <- asset.names
    er.gmin <- crossprod(w.gmin,er)
    sd.gmin <- sqrt(t(w.gmin) %*% cov.mat %*% w.gmin)
    gmin.port <- list("call" = call,
                      "er" = as.vector(er.gmin),
                      "sd" = as.vector(sd.gmin),
                      "weights" = w.gmin)
    class(gmin.port) <- "portfolio"
    gmin.port
  }

ir <- function(x) {
  mean(x)/sd(x)
}

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

calculate_portfolio_weights <- function(x, constrains_opt, nboot, factors = NULL) {
  w_estim <- optimalPortfolio(Sigma = cov(x), mu = colMeans(x), control = constrains_opt)
  w_bootparam <- michaud_parametric_bootstrap(x, B = nboot, option_list = constrains_opt) 
  w_boot <- michaud_bootstrap(x, B = nboot, option_list = constrains_opt) 
  k <- POETKhat(t(x))$K1HL
  w_factor_bootparam <- factor_parametric_bootstrap(x, B = nboot, n_factors  = k, option_list = constrains_opt, factors = NULL) 
  w_factor_boot <- factor_bootstrap(x, B = nboot, n_factors  = k, option_list = constrains_opt, factors = NULL) 
  if (is.list(factors)) {
    w_factor_bootparam_obsfactors1 <- factor_parametric_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors[[1]]) 
    w_factor_boot_obsfactors1 <- factor_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors[[1]]) 
    w_factor_bootparam_obsfactors2 <- factor_parametric_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors[[2]]) 
    w_factor_boot_obsfactors2 <- factor_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors[[2]]) 
    return(rbind(w_estim, w_bootparam, w_boot, w_factor_bootparam, w_factor_boot, w_factor_bootparam_obsfactors1, w_factor_boot_obsfactors1, w_factor_bootparam_obsfactors2, w_factor_boot_obsfactors2))
  } else{
    w_factor_bootparam_obsfactors <- factor_parametric_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors) 
    w_factor_boot_obsfactors <- factor_bootstrap(x, B = nboot, n_factors  = NULL, option_list = constrains_opt, factors) 
    return(rbind(w_estim, w_bootparam, w_boot, w_factor_bootparam, w_factor_boot, w_factor_bootparam_obsfactors, w_factor_boot_obsfactors))
  }
}