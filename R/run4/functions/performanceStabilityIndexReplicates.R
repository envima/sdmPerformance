#' @title Performance Stability Index Calculation
#'
#' @description
#' Computes a composite metric quantifying the performance and stability of spatial prediction models
#' using a combination of predictive accuracy (AUC, COR), predictive error measueres (MAE, BIAS), and inter-layer correlation (spatial stability)
#' across cross-validation folds.
#'
#' @param x A `terra::SpatRaster` object (multi-layer raster stack) representing prediction outputs across CV folds.
#' @param absence Optional. A `sf` object representing known absence locations. If not provided, `background` or `aa` must be used.
#' @param presence Required. A `sf` object representing known presence locations.
#' @param background Optional. Either a `sf` object with background points, or logical `TRUE` to auto-generate them. Logical `FALSE` to not calculate them
#' @param aa Optional. A `sf` object representing artificial absence points. If not provided, they are derived using AOA (Area of Applicability).
#' @param environmentalVariables Optional unless `background = TRUE` or `aa` is missing. A `terra::SpatRaster` with environmental covariates.
#' @param noPointsTesting Integer. Number of background or artificial absence points to generate.
#' @param prediction A `terra::SpatRaster` object with the prediction map.
#'
#' @return A named `list` with the following components:
#' \describe{
#'   \item{indexPA}{Performance index using provided absence data.}
#'   \item{indexPAA}{Performance index using artificially generated absence data (AOA-based).}
#'   \item{indexPBG}{Performance index using background data.}
#' }
#'
#' @details
#' The index combines accuracy metrics and a spatial correlation stability metric across folds to
#' estimate the reliability and generalizability of the model. COR and AUC measure predictive performance,
#' while layer-wise correlation (Kappa) measures spatial consistency across folds.
#'
#' @examples
#' \dontrun{
#'   result <- performanceStabilityIndex(x = prediction_stack,
#'                                       presence = presence_points,
#'                                       background = TRUE,
#'                                       environmentalVariables = env_rasters)
#'   result$indexPA
#' }
#'
#' @import terra sf dplyr CAST Metrics prg
#' @export

source("R/run4/functions/aacbs.R")
source("R/run3/functions/scale_metric.R")
source("R/run4/functions/sfbi.R")
source("R/run4/functions/confusionMatrix.R")

performanceStabilityIndex <- function(
    x = NA,
    prediction,
    presence = NA,
    absence = FALSE,
    background = TRUE,
    aa = TRUE,
    environmentalVariables = NA,
    noPointsTesting = NA,
    replicates=100
) {
  requireNamespace("sf")
  requireNamespace("terra")
  requireNamespace("CAST")
  requireNamespace("Metrics")
  requireNamespace("dplyr")
  
  if (!(inherits(x, "SpatRaster") || is.na(x)[1])) {stop("'x' must be a terra::SpatRaster object or NA.")}
  if (!inherits(prediction, "SpatRaster")) stop("'prediction' must be an spatRaster object.")
  if (!inherits(presence, "sf")) stop("'presence' must be an sf object.")
  if (!(inherits(absence, "sf") || isFALSE(absence))) stop("'absence' must be an sf object or FALSE")
  if (!(inherits(background, "sf") || is.logical(background))) stop("'background' must be an sf object or a logical (TRUE/FALSE)")
  if (!(inherits(aa, "sf") || is.logical(aa))) stop("'aa' must be an sf object or a logical (TRUE/FALSE)")
  if (!(inherits(environmentalVariables, "SpatRaster") || is.na(environmentalVariables))) stop("'environmentalVariables' must be either a terra::SpatRaster object or NA.")
  if (!(is.numeric(noPointsTesting) || is.na(noPointsTesting))) stop("'noPointsTesting' must be a numeric value or NA")
  
  
  
  if (inherits(absence, "sf") && nrow(absence) < 1) {
    absence <- FALSE
    message("Number of observations in absence is < 1. Set absence to FALSE. No metric on presence-absence data calculated.")
  }
  
  if (isFALSE(absence)[1] && isFALSE(background)[1] && isFALSE(aa)[1]) {
    stop("At least one of absence, background, or artificial absence (aa) must be provided.")
  }
  
  if ((isTRUE(background) ||  isTRUE(aa)) && !inherits(environmentalVariables, "SpatRaster")) {
    stop("Environmental variables must be provided to generate background or artificial absence data.")
  }
  
  if (is.logical(aa) && isTRUE(aa) && is.na(presence)[1]) {
    stop("Presence data must be provided to calculate artificial absence (AA) points.")
  }
  
  if (is.na(noPointsTesting)) {
    noPointsTesting <- nrow(presence)
  }
  
  #if ((isTRUE(aa) | isTRUE(background)) && noPointsTesting < 20){
  #  noPointsTesting<-20
  #  message("Number of artifical absence / background points to sample is set to 20.")
  #}
  
  #pred <- terra::mean(x)
  #pred=climateStability::rescale0to1(pred)
  
  if (inherits(x, "SpatRaster")) {
    stability <- stabilityRasters(x)
  } else {
    stability <- NA
    message("Skipping spatial stability calculation because 'x' is NA.")
  }
  
  
  
  
  
  if (isTRUE(background)) {
    message(paste("Start calculating metrics on presence-background data with ",replicates, "replicates."))
    indexPBG=do.call("rbind",lapply(1:replicates, function(i) {
      bg_df <- suppressMessages(as.data.frame(predicts::backgroundSample(environmentalVariables, n = noPointsTesting*5)))
      bg_df <- bg_df%>%dplyr::slice_sample( n=noPointsTesting)
      bg <- sf::st_as_sf(bg_df, coords = c("x", "y"), crs = terra::crs(environmentalVariables), remove = FALSE)
      
      inputPBG <- na.omit(rbind(
        data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
        data.frame(predicted = terra::extract(prediction, dplyr::sample_n(bg, noPointsTesting))[[2]], observed = 0)
      ))
      indexPBG <- indexCalculation(inputPBG, stability, prediction=prediction)
    })#end replicates
    )
    indexPBG=indexPBG %>% dplyr::summarize_all(mean, na.rm = TRUE)
    
    
  } else indexPBG <- NA
  
  if (!is.logical(absence) || !isFALSE(absence)) {
    inputPA <- na.omit(rbind(
      data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(prediction, absence)[[2]], observed = 0)
    ))
    indexPA <- indexCalculation(inputPA, stability, prediction=prediction)
  } else indexPA <- NA
  
  if (is.logical(aa) && isTRUE(aa)) {
    extr <- terra::extract(environmentalVariables, presence, ID = FALSE)
    aoa_result <- suppressMessages(CAST::aoa(newdata = environmentalVariables, train = extr, variables = "all",verbose = FALSE))
    # sample artificial absence points
    
    aa_mask <- aoa_result$AOA
    aa_mask[aa_mask > 0] <- NA
    # replicates
    message(paste("Start calculating metrics on presence-artificial-absence data with ",replicates, "replicates."))
    indexPAA=do.call("rbind",lapply(1:replicates, function(i) {
      aa_df <- suppressMessages(as.data.frame(predicts::backgroundSample(aa_mask, n = noPointsTesting*10, tryf = 30)))
      if(nrow(aa_df)>noPointsTesting){
        aa_df <- aa_df %>% dplyr::slice_sample(n = noPointsTesting)
      }
      aa <- sf::st_as_sf(aa_df, coords = c("x", "y"), crs = terra::crs(aa_mask), remove = FALSE)
      
      inputPAA <- na.omit(rbind(
        data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
        data.frame(predicted = terra::extract(prediction, dplyr::sample_n(aa, noPointsTesting))[[2]], observed = 0)
      ))
      indexPAA <- indexCalculation(inputPAA, stability, prediction=prediction)
      return(indexPAA)
    })#end replicates
    )
    indexPAA=indexPAA %>% dplyr::summarize_all(mean, na.rm = TRUE)
    
    rm(aa_mask, aoa_result,extr);gc()
  } else indexPAA <- NA
  
  
  if (!is.logical(absence) || !isFALSE(absence)) {
    inputPA <- na.omit(rbind(
      data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(prediction, absence)[[2]], observed = 0)
    ))
    indexPA <- indexCalculation(inputPA, stability, prediction=prediction)
  } else indexPA <- NA
  
  data=list(
    indexPA = indexPA,
    indexPAA = indexPAA,
    indexPBG = indexPBG
  )
  # Filter out NA data frames before binding
  data <- do.call(rbind, Filter(Negate(is.na), data))
  
  gc()
  return(data)
}

stabilityRasters <- function(r) {
  corr_matrix <- terra::layerCor(r, fun = "cor")$correlation
  corr_matrix[!lower.tri(corr_matrix, diag = FALSE)] <- NA
  valid_corrs <- na.omit(as.vector(corr_matrix))
  mean(valid_corrs)
}



indexCalculation <- function(inputDF, stability, prediction) {
  COR <- if (length(unique(inputDF$predicted)) > 1) cor(inputDF$observed, inputDF$predicted) else NA
  
  cm <- confusionMatrix(observation=inputDF$observed, predictions=inputDF$predicted)
  #cmPO <- confusionMatrix(observation=inputDF[inputDF$observed==1,]$observed, predictions=inputDF[inputDF$observed==1,]$predicted)
 # cmPO <- confusionMatrix_PO(observation=inputDF$observed, predictions=inputDF$predicted, prediction=prediction)
 # omissionRatePO<-omission(actual=cmPO$observed, predicted=cmPO$predicted)
  omissionRate<-omission(actual=cm$observed, predicted=cm$predicted)
  cmPBG <- confusionMatrix_PBG(observation=inputDF$observed, predictions=inputDF$predicted, prediction=prediction)
  omissionRatePBG<-omission(actual=cmPBG$observed, predicted=cmPBG$predicted)
  
  BS=brier_skill(observed=inputDF$observed, predicted = inputDF$predicted)
  BS_PO=brier_skill_po(observed=inputDF$observed, predicted = inputDF$predicted)
  FM= fowlkes_mallows(actual=cm$observed, predicted_binary=cm$predicted)
  MCC=mcc(actual=cm$observed, predicted_binary=cm$predicted)
  TopQ <- tcr(observed=inputDF$observed, predicted = inputDF$predicted, prediction=prediction)
  
  randomProbabilityValues<- terra::spatSample(prediction,size=5000,na.rm=T,as.df=TRUE)[[1]]
  boyce<-sfbi(prd1=inputDF[inputDF$observed==1,]$predicted, prd0=randomProbabilityValues, ktry=10)
  SBI_tp <-boyce[1]
  SBI_cr <-boyce[2]
  SBI_bs <-boyce[3]
  SBI_ps <-boyce[4]
  SBI_ad <-boyce[5]
  SBI_m <-boyce[6]
  
  FBP=Fbp(actual=cm$observed, predicted=cm$predicted)
  #FBP=Fbp(actual=inputDF$observed, predicted=inputDF$predicted)
  SEDI= sedi(actual=cm$observed, predicted=cm$predicted)
  ORSS=orss(actual=cm$observed, predicted=cm$predicted)
  AACBS=AACBS(predicted=inputDF$predicted, observed=inputDF$observed)
  
  
  AUC <- Metrics::auc(inputDF$observed, inputDF$predicted)
  PRG <- prg::calc_auprg(prg::create_prg_curve(inputDF$observed, inputDF$predicted))
  MAE <- Metrics::mae(inputDF$observed, inputDF$predicted)
  BIAS <-  abs(Metrics::bias(inputDF$observed, inputDF$predicted))
  evalDat=mecofun::evalSDM(inputDF$observed, inputDF$predicted)
  TSS=evalDat$TSS
  Kappa=evalDat$Kappa
  PCC=evalDat$PCC
  Sens=evalDat$Sens
  Spec=evalDat$Spec
 
  
  if(is.na(stability)){
    metric= mean(c(Spec, COR,PCC,1-MAE), na.rm=T)
  } else {
    metric= mean(c(Spec, COR,PCC,1-MAE,stability),na.rm=T)
  }
  
  # all metrics in one df
  result=data.frame(metric=metric,Fbp=FBP,AACBS=AACBS,BS=BS,omissionRatePBG=omissionRatePBG,omissionRate=omissionRate,SBI_tp=SBI_tp, SBI_cr=SBI_cr, SBI_bs=SBI_bs, SBI_ps=SBI_ps, SBI_ad=SBI_ad, SBI_m=SBI_m, TopQ=TopQ,BS_PO=BS_PO,FM=FM,MCC=MCC, SEDI=SEDI,ORSS=ORSS,AUC=AUC, COR=COR, Spec,Sens,Kappa,PCC, TSS,stability=stability, PRG=PRG, MAE=MAE, BIAS=BIAS, noPresencePoints=nrow(inputDF[inputDF$observed ==1,]))
  return(result)
}



#modiefied from Metrics::f1
Fbp <- function (actual, predicted) {
  act <- unique(actual)
  pred <- unique(predicted)
  
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  tn <- sum(actual == 0 & predicted == 0)
  
  if (tp == 0) {
    return(0)
  } else {
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    #precision2 = precision/(1-precision)
    fbp=(2*tp)/(tp+fn+fp)
    return(fbp)
  }
}


# Symmetric Extremal Dependence Index (SEDI)
# eps is a very small number used to prevent division by zero or taking the logarithm of 0 or 1.
sedi <- function (actual, predicted, eps = 1e-10) {
  act <- unique(actual)
  pred <- unique(predicted)
  
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  tn <- sum(actual == 0 & predicted == 0)
  
  
  #precision <- tp / (tp + fp)
  H <- tp / (tp + fn)
  FPR=fp/(fp+tn)
  
  # Guard against 0 or 1 (SEDI undefined at boundaries)
  H <- pmin(pmax(H, eps), 1 - eps)
  FPR <- pmin(pmax(FPR, eps), 1 - eps)
  
  sedi_val <- (log(FPR) - log(H) + log(1 - FPR) - log(1 - H)) /
    (log(FPR) + log(H) + log(1 - H) + log(1 - FPR))
  
  return(sedi_val)
  
}


#Odds Ratio Skill Score (ORSS; Stephenson 2000)
# eps is a very small number used to prevent division by zero or taking the logarithm of 0 or 1.
orss <- function (actual, predicted) {
  act <- unique(actual)
  pred <- unique(predicted)
  
  tp <- sum(actual == 1 & predicted == 1)
  fp <- sum(actual == 0 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  tn <- sum(actual == 0 & predicted == 0)
  
  
  a=tp
  b=fp
  c=fn
  d=tn
  
  ORSS = (a*d - b*c)/(a*d + b*c)
  return(ORSS)
}


omission <- function (actual, predicted) {
  act <- unique(actual)
  pred <- unique(predicted)
  
  tp <- sum(actual == 1 & predicted == 1)
  fn <- sum(actual == 1 & predicted == 0)
  
  
  omissionRate=fn/(tp+fn)
  
  return(omissionRate)
}

brier_score <- function(actual, predicted_probs) {
  if(length(actual) != length(predicted_probs)) stop("Vectors must be same length")
  mean((predicted_probs - actual)^2)
}

#Fowlkes-Mallows Index
fowlkes_mallows <- function(actual, predicted_binary) {
  tp <- sum(actual == 1 & predicted_binary == 1)
  fp <- sum(actual == 0 & predicted_binary == 1)
  fn <- sum(actual == 1 & predicted_binary == 0)
  
  if(tp == 0) return(0)
  sqrt((tp / (tp + fp)) * (tp / (tp + fn)))
}

# Matthews Correlation Coefficient (MCC)
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



#Top-q% capture rate
#tcr <- function(predicted, observed, q = 0.10) {
#  stopifnot(length(predicted) == length(observed))
#  pres <- predicted[observed == 1]
#  bg   <- predicted[observed == 0]
#  thr  <- quantile(bg, probs = 1 - q, type = 7)
#  mean(pres >= thr)
#}





#Top-q% capture rate
tcr <- function(predicted, observed, prediction,q = 0.10) {
  stopifnot(length(predicted) == length(observed))
  pres <- predicted[observed == 1]
  prediction_values=terra::values(prediction,na.rm=T)
  thr  <- quantile( prediction_values, probs = 1 - q, type = 7)
  mean(pres >= thr)
}
