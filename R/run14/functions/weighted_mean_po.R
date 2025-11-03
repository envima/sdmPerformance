# Weighted mean function
weighted_extreme_mean <- function(TopQ, omissionRate, alpha = 2) {
  # values: numeric vector of metrics (between 0 and 1)
  # alpha: controls how much extreme values are emphasized
  #weight=abs(TopQ-omissionRate)
  values=c(TopQ, omissionRate)
  # min=TopQ
  #  max=omissionRate
  
  
  # weighted_mean=TopQ*weight + omissionRate*1-weight
  
  # if (weighted_mean > omissionRate) weighted_mean <-omissionRate
  # if (weighted_mean < TopQ) weighted_mean <-TopQ
  
  
  #TopQ
  #weighted_mean
  #omissionRate
  
  
  weights <- 1 + alpha * abs(values - 0.5)
  
  weighted_mean <- sum(values * weights) / sum(weights)
  
  #if (omissionRate < 0.5) weighted_mean <- mean(c(omissionRate,weighted_mean))
  
  # Step 2: Compute dynamic penalty for low omissionRate
  # Low omissionRate -> strong penalty, high omissionRate -> negligible penalty
  penalty_factor <- 1 - (1 - omissionRate)^2
  
  # Step 3: Apply penalty
  weighted_mean <- weighted_mean * penalty_factor
  
  
  return(weighted_mean)
}
