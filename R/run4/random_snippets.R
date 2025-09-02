
testDF=data.frame(predicted=c(predicted[idx],runif(n = length(predicted[idx]),min=0,max=0.5)), observed=c(runif(n = length(predicted[idx]),min=1,max=1), runif(n = length(predicted[idx]),min=0,max=0)))
testDF=data.frame(predicted=c(predicted[idx],runif(n = length(predicted[idx]),min=0,max=0.5)), observed=c(runif(n = length(predicted[idx]),min=1,max=1), runif(n = length(predicted[idx]),min=0,max=0.5)))



brier_skill(observed=testDF$observed, predicted=testDF$predicted)
brier_skill(observed=inputDF$observed, predicted = inputDF$predicted)

brier_skill_po(observed=inputDF$observed, predicted = inputDF$predicted)



testi=data %>% dplyr::filter(species=="VS01", points==120, testData==1, size=="KNNDM", replicate==1)




AACBS(observed=testDF$observed, predicted=testDF$predicted)
AACBS(observed=inputDF$observed, predicted = inputDF$predicted)



presence_rank_score <- function(predicted, observed, spatial_weight = NULL) {
  stopifnot(length(predicted) == length(observed))
  idx <- observed == 1
  if (!any(idx)) stop("No presences in 'observed'.")
  
  ranks <- rank(predicted[idx]) / sum(idx)
  if (!is.null(spatial_weight)) ranks <- ranks * spatial_weight[idx]
  
  score <- 1 - mean((ranks - 1)^2)
  return(score)
}



# Only consider presence points
idx <- observed == 1
pred_pres <- predicted[idx]
n <- length(pred_pres)

# --- Rank component ---
ranks <- rank(pred_pres) / n    # normalize rank to [0,1]

# --- Scaled probability component ---
prob_scaled <- (pred_pres - min(pred_pres)) / (max(pred_pres) - min(pred_pres))

# --- Combine both components equally ---
hybrid_score <- mean(ranks)



# Top Q -> minimum of the metric
TopQ=data$TopQ[i]
omissionRate=data$omissionRate[i]
# Max 1-omission rate, calculate the threshold here with fbp to aviod absence metrics
omissionRate=1-omissionRate
#this is equivalent to sensititvity

# most likely rnage of the metric:
r=c(TopQ, omissionRate)
mean(c(TopQ, omissionRate))

pcc=data$PCC_scaled[i]
cor=data$COR_scaled[i]
mean(c(pcc,cor))


sens of 0.75 means between 0 and 0.75
top q of 0.75 means between 0.75 and 1

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
  
  # Step 4: Dynamic boost for high TopQ
  #boost_factor <- (1 - (1 - TopQ)^0.5)
  
  
  #weighted_mean_Boost <- weighted_mean * (1-boost_factor*0.5) + TopQ*boost_factor*0.5
  #weighted_mean <- weighted_mean  + TopQ*boost_factor*0.2
  
  #boost <- (TopQ^2 - 0.5^2) #* 0.1  # small additive boost, only significant if TopQ is high
  #weighted_mean <- weighted_mean + boost
  
  return(weighted_mean)
}

# Example usage
sens <- omissionRate
topq_sens <- TopQ

metrics <- c(sens, topq_sens)
result <- weighted_extreme_mean(metrics, alpha = 2)
print(result)
