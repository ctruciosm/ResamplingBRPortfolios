#############################
### Performance Measures  ###
#############################


medidas <- function(x, rf = 0.5) {
  # Annualized Average
  AV <- mean(x)
  # Annualized SD
  SD <- sd(x)
  # Information (or Sharpe) Ratio
  #IR <- AV/SD
  SR <- (mean(x) - rf)/sd(x)
  # Adjusted Sharpe Ratio
  ASR <- SR*(1 + (moments::skewness(x)/6)*SR - ((moments::kurtosis(x) - 3)/24)*SR^2)
  # Sortino Ratio
  #SO <- AV/sqrt(mean(ifelse(x < 0, 0, x^2)))
  SO <- (mean(x) - rf)/sqrt(mean(ifelse(x - rf < 0, 0, (x - rf)^2)))
  output <- c(12*AV, sqrt(12)*SD, SR, ASR, SO)
  return(output)
}


calculate_to_sspw <- function(previous_weights, desired_weights, oos_returns, p) {
  l <- length(previous_weights)/p
  to <- rep(NA, l)
  sspw <- rep(NA, l)
  for (k in 1:l) {
    num <- previous_weights[(p*(k - 1) + 1):(p*k)]*(1 + oos_returns/100)
    den <- sum(num)
    updated_weights <- num/den
    to[k] <- sum(abs(desired_weights[(p*(k - 1) + 1):(p*k)] - updated_weights))
    sspw[k] <- sum(desired_weights[(p*(k - 1) + 1):(p*k)]^2)
  }
  return(list(to, sspw))
}