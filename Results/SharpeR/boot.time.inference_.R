boot.time.inference_ <- function (ret, b, M=4999, Delta.null = 0, digits = 4) 
{
    T = length(ret[, 1])
    l = floor(T/b)
    Delta.hat = sharpe.ratio.diff(ret)
    d = abs(Delta.hat - Delta.null)/compute.se.Parzen.pw(ret)

    library(parallel)
    library(MASS)
    trials <- seq(1, M)

    fx <- function(trial){
        pvalue = 0
        ret.star = ret[cbb.sequence(T, b), ]
        Delta.hat.star = sharpe.ratio.diff(ret.star)
        ret1.star = ret.star[, 1]
        ret2.star = ret.star[, 2]
        mu1.hat.star = mean(ret1.star)
        mu2.hat.star = mean(ret2.star)
        gamma1.hat.star = mean(ret1.star^2)
        gamma2.hat.star = mean(ret2.star^2)
        gradient = rep(0, 4)
        gradient[1] = gamma1.hat.star/(gamma1.hat.star - mu1.hat.star^2)^1.5
        gradient[2] = -gamma2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5
        gradient[3] = -0.5 * mu1.hat.star/(gamma1.hat.star - mu1.hat.star^2)^1.5
        gradient[4] = 0.5 * mu2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5
        y.star = data.frame(ret1.star - mu1.hat.star, ret2.star - 
            mu2.hat.star, ret1.star^2 - gamma1.hat.star, ret2.star^2 - gamma2.hat.star)
        Psi.hat.star = matrix(0, 4, 4)
        for (j in (1:l)) {
            zeta.star = b^0.5 * colMeans(y.star[((j - 1) * b + 1):(j * b), ])
            Psi.hat.star = Psi.hat.star + zeta.star %*% t(zeta.star)
        }
        Psi.hat.star = Psi.hat.star/l
        se.star = as.numeric(sqrt(t(gradient) %*% Psi.hat.star %*% gradient/T))
        d.star = abs(Delta.hat.star - Delta.hat)/se.star
        if (d.star >= d) { pvalue = 1 }
        return(pvalue)
    }

    numCores <- detectCores()
    results <- mclapply(trials, fx, mc.cores = numCores)
    p.value = ( sum(unlist(results)) + 1 ) / (M + 1)
    list( Difference = round(Delta.hat, digits), p.Value = round(p.value, digits) )
}

