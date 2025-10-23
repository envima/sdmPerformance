

# Presence-only Brier Skill Score (baseline = 0.5 at presence sites)
brier_skill <- function(observed, predicted, clip01 = TRUE) {
  stopifnot(length(observed) == length(predicted))
 # if (!all(observed %in% c(0,1))) stop("'observed' must be 0/1.")
  if (any(predicted < 0 | predicted > 1)) stop("'predicted' must be in [0,1].")
  idx  <- observed == 1
  if (!any(idx)) stop("No presences in 'observed'.")
  bs_po <- mean((predicted[idx] - 1)^2)
  brier_score_all <- mean((predicted - observed)^2)
 
  max_val=0.25 # ~ assumed random prediction brier score  # 0.25 = (0.5 - 1)^2
  min_val=0
  value=brier_score_all
#  if(value < 0.1) value=0.1 
  scaled <- 1-(value - min_val) / (max_val - min_val)
 # scaled <- (value - min_val) / (max_val - min_val)
  
  
  #bss  <- 1 - bs_po / 0.25                   
  #if (clip01) bss <- max(0, min(1, bss))
  return(scaled)
}




# Presence-only Brier Skill Score (baseline = 0.5 at presence sites)
brier_skill_po <- function(observed, predicted, clip01 = TRUE) {
  stopifnot(length(observed) == length(predicted))
 # if (!all(observed %in% c(0,1))) stop("'observed' must be 0/1.")
  if (any(predicted < 0 | predicted > 1)) stop("'predicted' must be in [0,1].")
  idx  <- observed == 1
  if (!any(idx)) stop("No presences in 'observed'.")
# df_bs=data.frame(predicted=c(predicted[idx],runif(n = length(predicted[idx]),min=0,max=0)), observed=c(runif(n = length(predicted[idx]),min=1,max=1), runif(n = length(predicted[idx]),min=0,max=0.5)))
  
  
  
  #bs_po <- mean((df_bs$predicted - df_bs$observed)^2)
  bs_po <- mean((predicted[idx] - 1)^2)
  
  
  max_val=0.25 # ~ assumed random prediction brier score  # 0.25 = (0.5 - 1)^2
  min_val=0
  value=bs_po
 # if(value < 0.01) value=0.01 
  scaled <- 1-(value - min_val) / (max_val - min_val)
  
  #scaled<-bs_po
  #bss  <- 1 - bs_po / 0.25                   
  if (clip01) scaled <- max(0, min(1, scaled))
  return(scaled)
}

# artifical absence corrected breier score
# Updated composite with adaptive use of presence-only Brier skill and PB-based TSS/COR
AACBS <- function(observed, predicted, optimize_tss = TRUE, thresholds = NULL) {
  stopifnot(length(observed) == length(predicted))
 # if (!all(observed %in% c(0,1))) stop("'observed' must be 0/1.")
  if (any(predicted < 0 | predicted > 1)) stop("'predicted' must be in [0,1].")
  
  # Correlation (scaled via your custom function)
  cor_val <- suppressWarnings(cor(predicted, observed, method = "pearson"))
  cor_val <- scale_metric(cor_val, "COR")
  cor_val <- ifelse(is.na(cor_val), 0, cor_val)
  
  # TSS on presence–background
  if (isTRUE(optimize_tss)) {
    thr <- if (is.null(thresholds)) sort(unique(predicted)) else sort(unique(thresholds))
    best_tss <- -Inf
    for (t in thr) {
      pred_bin <- as.integer(predicted >= t)
      TP <- sum(pred_bin == 1 & observed == 1)
      TN <- sum(pred_bin == 0 & observed == 0)
      FP <- sum(pred_bin == 1 & observed == 0)
      FN <- sum(pred_bin == 0 & observed == 1)
      sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
      tss  <- sens + spec - 1
      if (tss > best_tss) best_tss <- tss
    }
    tss_val <- best_tss
  } else {
    pred_bin <- as.integer(predicted >= 0.5)
    TP <- sum(pred_bin == 1 & observed == 1)
    TN <- sum(pred_bin == 0 & observed == 0)
    FP <- sum(pred_bin == 1 & observed == 0)
    FN <- sum(pred_bin == 0 & observed == 1)
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
    tss_val <- sens + spec - 1
  }
  tss_val <- scale_metric(tss_val, "TSS")     # map to [0,1] as in your setup
  
  # Brier skill scores
  #bs <- brier_skill_po(observed, predicted, clip01 = TRUE) # presence-only skill
  bs <- brier_skill(observed, predicted, clip01 = TRUE) 
  
  
  # Hybrid correction: if presence-only skill is high, add TSS and COR
  # weight_factor grows as presence-only BS improves (i.e., as bs increases)
  # e.g., with bs = 0.64 (equiv. BS_pres = 0.09), then weight_factor = 0.64
  #weight_factor <- bs*0.5
  #score    <-  (weight_factor)*bs + (tss_val + cor_val)*0.5*(1-weight_factor)
  score=mean(c(bs,cor_val,tss_val))
  
  # Base is presence-only Brier skill; add correction against over-optimism
  #score <- bs + correction
  
  # Optional clipping for comparability
  #score <- max(0, min(1, score))
  as.numeric(score)
  return(score)
}


