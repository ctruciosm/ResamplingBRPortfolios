block.size.calibrate_ <- function (ret, b.vec = c(1L, 3L, 6L, 10L), alpha = 0.05, M = 199, K = 1000, b.av = 5, T.start = 50) {
    
    b.len = length(b.vec)
    Delta.hat = sharpe.ratio.diff(ret)
    ret1 = ret[, 1]
    ret2 = ret[, 2]
    T = length(ret1)
    
    fit1 = lm(ret1[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    fit2 = lm(ret2[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
    coef1 = as.numeric(fit1$coef)
    coef2 = as.numeric(fit2$coef)
    resid.mat = cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))
    
    library(parallel)
    library(MASS)
    trials <- seq(1, K)
    
    fx <- function(trial){     
        emp.reject.probs = rep(0, b.len)
        Var.data = matrix(0, T.start + T, 2)
        Var.data[1, 1] = ret[1, 1]
        Var.data[1, 2] = ret[1, 2]
        
        resid.mat.star = rbind(c(0, 0), resid.mat[sb.sequence(T - 1, b.av, T.start + T - 1), ])
        for (t in (2:(T.start + T))) {
            Var.data[t, 1] = coef1[1] + coef1[2] * Var.data[t - 1, 1] + coef1[3] * Var.data[t - 1, 2] + resid.mat.star[t, 1]
            Var.data[t, 2] = coef2[1] + coef2[2] * Var.data[t - 1, 1] + coef2[3] * Var.data[t - 1, 2] + resid.mat.star[t, 2]
        }
        Var.data.trunc = Var.data[(T.start + 1):(T.start + T), ]
        for (j in (1:b.len)) {
            p.Value = boot.time.inference(Var.data.trunc, b.vec[j], M, Delta.hat)$p.Value          
            if (p.Value <= alpha) { emp.reject.probs[j] = emp.reject.probs[j] + 1 }
        }
        return(emp.reject.probs)
    }

    numCores <- detectCores()
    results <- mclapply(trials, fx, mc.cores = numCores)
    emprejectprobslist = t(array(unlist(results), dim=c(b.len, K)))
    emprejectprobs = colSums( emprejectprobslist ) / K
    b.order = order(abs(emprejectprobs - alpha))
    b.opt = b.vec[b.order[1]]
    b.vec.with.probs = rbind(b.vec, emprejectprobs)
    colnames(b.vec.with.probs) = rep("", length(b.vec))
    list(Empirical.Rejection.Probs = b.vec.with.probs, b.optimal = b.opt)
}

