#' @title Calculate SDM Metrics
#'
#' @description
#' Computes a comprehensive set of species distribution model (SDM) evaluation metrics
#' including correlation, AUC, PRG, MAE, BIAS, and several skill scores.
#'
#' @param inputDF A data frame with columns \code{predicted} and \code{observed}.
#' @param prediction A \code{terra::SpatRaster} object containing predicted values.
#'
#' @return A \code{data.frame} containing multiple SDM evaluation metrics for the input data.
#' @export


source("R/run8/functions/sfbi.R")
source("R/run8/functions/confusionMatrix.R")


indexCalculation <- function(inputDF, prediction) {
  
  # --------------------------------------------------------------------------
  # Compute correlation between observed and predicted
  # --------------------------------------------------------------------------
  COR <- if (length(unique(inputDF$predicted)) > 1) cor(inputDF$observed, inputDF$predicted) else NA
  
  # --------------------------------------------------------------------------
  # Generate confusion matrix for observed vs predicted
  # --------------------------------------------------------------------------
  cm <- confusionMatrix(observation = inputDF$observed, predictions = inputDF$predicted)
  
  # Compute omission rate from confusion matrix
  omissionRate <- omission(actual = cm$observed, predicted = cm$predicted)
  
  # --------------------------------------------------------------------------
  # Sample random predictions from raster for Boyce index / stability metrics
  # --------------------------------------------------------------------------
  randomProbabilityValues <- terra::spatSample(prediction, size = 5000, na.rm = TRUE, as.df = TRUE)[[1]]
  
  # Compute Boyce index for presence vs random predictions
  boyce <- sfbi(prd1 = inputDF[inputDF$observed == 1,]$predicted,
                prd0 = randomProbabilityValues, ktry = 10)
  SBI_tp <- boyce[1]; SBI_cr <- boyce[2]; SBI_bs <- boyce[3]
  SBI_ps <- boyce[4]; SBI_ad <- boyce[5]; SBI_m <- boyce[6]
  
  # --------------------------------------------------------------------------
  # Compute other skill metrics
  # --------------------------------------------------------------------------
  FBP <- Fbp(actual = cm$observed, predicted = cm$predicted)       # Modified F1 score
  SEDI <- sedi(actual = cm$observed, predicted = cm$predicted)    # Symmetric Extremal Dependence Index
  ORSS <- orss(actual = cm$observed, predicted = cm$predicted)    # Odds Ratio Skill Score
  
  # --------------------------------------------------------------------------
  # Standard accuracy metrics
  # --------------------------------------------------------------------------
  AUC <- Metrics::auc(inputDF$observed, inputDF$predicted)        # Area Under Curve
  PRG <- prg::calc_auprg(prg::create_prg_curve(inputDF$observed, inputDF$predicted)) # Precision-Recall Gain
  MAE <- Metrics::mae(inputDF$observed, inputDF$predicted)        # Mean Absolute Error
  BIAS <- abs(Metrics::bias(inputDF$observed, inputDF$predicted)) # Absolute bias
  
  # --------------------------------------------------------------------------
  # Calculate additional evaluation metrics using mecofun
  # --------------------------------------------------------------------------
  evalDat <- mecofun::evalSDM(inputDF$observed, inputDF$predicted)
  TSS <- evalDat$TSS; Kappa <- evalDat$Kappa; PCC <- evalDat$PCC
  Sens <- evalDat$Sens; Spec <- evalDat$Spec
  
  # --------------------------------------------------------------------------
  # Combine all metrics into a single data frame
  # --------------------------------------------------------------------------
  result <- data.frame(
    Fbp = FBP, omissionRate = omissionRate,
    SBI_tp = SBI_tp, SBI_cr = SBI_cr, SBI_bs = SBI_bs,
    SBI_ps = SBI_ps, SBI_ad = SBI_ad, SBI_m = SBI_m,
    SEDI = SEDI, ORSS = ORSS, AUC = AUC, COR = COR,
    Spec = Spec, Sens = Sens, Kappa = Kappa, PCC = PCC, TSS = TSS,
    PRG = PRG, MAE = MAE, BIAS = BIAS,
    noPresencePoints = nrow(inputDF[inputDF$observed == 1,])
  )
  
  return(result)
}





# ============================================================================
# Metric helper functions
# ============================================================================

#' @title Fbp: Modified F1 Score
#' @description Computes Fbp (Fowlkes-F1) score from binary predictions.
#' @param actual Vector of observed binary outcomes.
#' @param predicted Vector of predicted binary outcomes.
#' @return Numeric Fbp score.
#' @export
Fbp <- function(actual, predicted) {
  # Compute true positives, false positives, false negatives
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  
  if (tp == 0) return(0)
  
  # Compute modified F1 / Fbp
  fbp <- (2 * tp) / (tp + fn + fp)
  return(fbp)
}

#' @title Symmetric Extremal Dependence Index (SEDI)
#' @description Evaluates prediction skill for rare events.
#' @param actual Vector of observed binary outcomes.
#' @param predicted Vector of predicted binary outcomes.
#' @param eps Small number to prevent division by zero. Default = 1e-10.
#' @return Numeric SEDI value.
#' @export
sedi <- function(actual, predicted, eps = 1e-10) {
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  tn <- sum(actual == 0 & predicted == 0)
  
  # Calculate hit rate and false positive rate
  H <- pmin(pmax(tp / (tp + fn), eps), 1 - eps)
  FPR <- pmin(pmax(fp / (fp + tn), eps), 1 - eps)
  
  # Compute SEDI
  sedi_val <- (log(FPR) - log(H) + log(1 - FPR) - log(1 - H)) /
    (log(FPR) + log(H) + log(1 - H) + log(1 - FPR))
  return(sedi_val)
}

#' @title Odds Ratio Skill Score (ORSS)
#' @description Evaluates prediction skill based on contingency table odds ratio.
#' @param actual Vector of observed binary outcomes.
#' @param predicted Vector of predicted binary outcomes.
#' @return Numeric ORSS value.
#' @export
orss <- function(actual, predicted) {
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  tn <- sum(actual == 0 & predicted == 0)
  
  ORSS <- (tp * tn - fp * fn) / (tp * tn + fp * fn)
  return(ORSS)
}

#' @title Omission Rate
#' @description Computes omission rate for binary predictions.
#' @param actual Vector of observed binary outcomes.
#' @param predicted Vector of predicted binary outcomes.
#' @return Numeric omission rate.
#' @export
omission <- function(actual, predicted) {
  tp <- sum(actual == 1 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  return(fn / (tp + fn))
}

#' @title Brier Score
#' @description Measures mean squared error of probabilistic predictions.
#' @param actual Vector of observed outcomes (0/1).
#' @param predicted_probs Vector of predicted probabilities.
#' @return Numeric Brier score.
#' @export
brier_score <- function(actual, predicted_probs) {
  if(length(actual) != length(predicted_probs)) stop("Vectors must be same length")
  mean((predicted_probs - actual)^2)
}

#' @title Fowlkes-Mallows Index
#' @description Computes Fowlkes-Mallows index for binary predictions.
#' @param actual Vector of observed binary outcomes.
#' @param predicted_binary Vector of predicted binary outcomes.
#' @return Numeric FM index.
#' @export
fowlkes_mallows <- function(actual, predicted_binary) {
  tp <- sum(actual == 1 & predicted_binary == 1)
  fp <- sum(actual == 0 & predicted_binary == 1)
  fn <- sum(actual == 1 & predicted_binary == 0)
  
  if(tp == 0) return(0)
  sqrt((tp / (tp + fp)) * (tp / (tp + fn)))
}

#' @title Matthews Correlation Coefficient (MCC)
#' @description Computes MCC for binary classification.
#' @param actual Vector of observed binary outcomes.
#' @param predicted_binary Vector of predicted binary outcomes.
#' @return Numeric MCC value.
#' @export
mcc <- function(actual, predicted_binary) {
  tp <- sum(actual == 1 & predicted_binary == 1)
  tn <- sum(actual == 0 & predicted_binary == 0)
  fp <- sum(actual == 0 & predicted_binary == 1)
  fn <- sum(actual == 1 & predicted_binary == 0)
  
  numerator <- tp * tn - fp * fn
  denominator <- sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
  if(denominator == 0) return(0)
  numerator / denominator
}

#' @title Top-q% Capture Rate
#' @description Computes proportion of presences captured in top q% of predicted values.
#' @param predicted Vector of predicted values.
#' @param observed Vector of observed binary outcomes.
#' @param prediction \code{terra::SpatRaster} object used to compute threshold.
#' @param q Numeric. Top proportion threshold (default 0.10 = top 10%).
#' @return Numeric proportion of presences captured.
#' @export
tcr <- function(predicted, observed, prediction, q = 0.10) {
  stopifnot(length(predicted) == length(observed))
  pres <- predicted[observed == 1]
  
  # Compute threshold value from raster
  prediction_values <- terra::values(prediction, na.rm = TRUE)
  thr <- quantile(prediction_values, probs = 1 - q, type = 7)
  
  # Return proportion of presences above threshold
  mean(pres >= thr)
}