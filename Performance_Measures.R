#############################
### Performance Measures  ###
#############################


medidas <- function(x, MAR = 0) {
  # Annualized Average
  AV <- mean(x)
  # Annualized SD
  SD <- sd(x)
  # Information (or Sharpe) Ratio
  #IR <- AV/SD
  SR <- mean(x)/sd(x)
  # Adjusted Sharpe Ratio
  ASR <- SR*(1 + (moments::skewness(x)/6)*SR - ((moments::kurtosis(x) - 3)/24)*SR^2)
  # Sortino Ratio
  #SO <- AV/sqrt(mean(ifelse(x < 0, 0, x^2)))
  SO <- mean(x)/sqrt(mean(ifelse(x < 0, 0, x^2)))
  output <- c(12*AV, sqrt(12)*SD, SR, ASR, SO)
  return(output)
}

weights_measures <- function(previous_weights, desired_weights, oos_returns) {
  dim_aux <- ncol(desired_weights)
  if (is.null(dim_aux)) {
    num <- previous_weights*(1 + oos_returns/100)
    den <- sum(num)
    updated_weights <- num/den  
    to <- sum(abs(desired_weights - updated_weights))
    sspw <- sum(desired_weights^2)
  } else{
    to = rep(0, dim_aux)
    sspw = rep(0, dim_aux)
    for (i in 1:dim_aux) {
      num <- previous_weights[,i]*(1 + oos_returns/100)
      den <- sum(num)
      updated_weights <- num/den  
      to[i] <- sum(abs(desired_weights[,i] - updated_weights))
      sspw[i] <- sum(desired_weights[,i]^2)
    }
  }
  return(list(to, sspw))
}