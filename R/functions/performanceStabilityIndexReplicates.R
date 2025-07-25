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
#' @param trainingData Presence only data on which the model was trained. A `sf` object.
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

performanceStabilityIndex <- function(
    x = NA,
    prediction,
    presence = NA,
    absence = FALSE,
    background = TRUE,
    aa = TRUE,
    environmentalVariables = NA,
    trainingData = NA,
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
  if (!(inherits(trainingData, "sf") || is.na(trainingData))) stop("'trainingData' must be an sf object or NA")
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
  
  if (is.logical(aa) && isTRUE(aa) && is.na(trainingData)[1]) {
    stop("Training data must be provided to calculate artificial absence (AA) points.")
  }
  
  if (is.na(noPointsTesting)) {
    noPointsTesting <- nrow(presence)
  }
  
  if ((isTRUE(aa) | isTRUE(background)) && noPointsTesting < 20){
    noPointsTesting<-20
    message("Number of artifical absence / background points to sample is set to 20.")
  }
  
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
      indexPBG <- indexCalculation(inputPBG, stability)
    })#end replicates
    )
    indexPBG=indexPBG %>% dplyr::summarize_all(mean, na.rm = TRUE)
    
    
  } else indexPBG <- NA
  
  if (!is.logical(absence) || !isFALSE(absence)) {
    inputPA <- na.omit(rbind(
      data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(prediction, absence)[[2]], observed = 0)
    ))
    indexPA <- indexCalculation(inputPA, stability)
  } else indexPA <- NA
  
  if (is.logical(aa) && isTRUE(aa)) {
    extr <- terra::extract(environmentalVariables, trainingData, ID = FALSE)
    aoa_result <- suppressMessages(CAST::aoa(newdata = environmentalVariables, train = extr, variables = "all",verbose = FALSE))
    # sample artificial absence points
    
    aa_mask <- aoa_result$AOA
    aa_mask[aa_mask > 0] <- NA
    # replicates
    message(paste("Start calculating metrics on presence-artificial-absence data with ",replicates, "replicates."))
    indexPAA=do.call("rbind",lapply(1:replicates, function(i) {
      aa_df <- suppressMessages(as.data.frame(predicts::backgroundSample(aa_mask, n = noPointsTesting*10, tryf = 10)))
      if(nrow(aa_df)>noPointsTesting){
        aa_df <- aa_df %>% dplyr::slice_sample(n = noPointsTesting)
      }
      aa <- sf::st_as_sf(aa_df, coords = c("x", "y"), crs = terra::crs(aa_mask), remove = FALSE)
      
      inputPAA <- na.omit(rbind(
        data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
        data.frame(predicted = terra::extract(prediction, dplyr::sample_n(aa, noPointsTesting))[[2]], observed = 0)
      ))
      indexPAA <- indexCalculation(inputPAA, stability)
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
    indexPA <- indexCalculation(inputPA, stability)
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



#indexCalculation <- function(inputDF, stability) {
#  COR <- if (length(unique(inputDF$predicted)) > 1) cor(inputDF$observed, inputDF$predicted) else NA
#  if (is.na(COR)) {
#    COR <- NA
#  } else if (COR < 0) {
#    COR <- 0
#  }
#  AUC <- Metrics::auc(inputDF$observed, inputDF$predicted)
#  PRG <- prg::calc_auprg(prg::create_prg_curve(inputDF$observed, inputDF$predicted))
#  MAE <- Metrics::mae(inputDF$observed, inputDF$predicted)
#  BIAS <-  abs(Metrics::bias(inputDF$observed, inputDF$predicted))
#  evalDat=mecofun::evalSDM(inputDF$observed, inputDF$predicted)
#  TSS=evalDat$TSS
#  Kappa=evalDat$Kappa
#  PCC=evalDat$PCC
#  Sens=evalDat$Sens
#  Spec=evalDat$Spec
#  
#  if(is.na(stability)){
#    metric= mean(c(COR,AUC,1-MAE,1-BIAS), na.rm=T)
#  } else {
#    metric= mean(c(COR,stability,AUC,1-MAE,1-BIAS),na.rm=T)
#  }
#  
#  # all metrics in one df
#  result=data.frame(metric=metric, AUC=AUC, COR=COR, Spec,Sens,Kappa,PCC, TSS,stability=stability, PRG=PRG, MAE=MAE, BIAS=BIAS, noPresencePoints=nrow(inputDF[inputDF$observed ==1,]))
#  return(result)
#}



indexCalculation <- function(inputDF, stability) {
  COR <- if (length(unique(inputDF$predicted)) > 1) cor(inputDF$observed, inputDF$predicted) else NA
  if (is.na(COR)) {
    COR <- NA
  } else if (COR < 0) {
    COR <- 0
  }
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
  result=data.frame(metric=metric, AUC=AUC, COR=COR, Spec,Sens,Kappa,PCC, TSS,stability=stability, PRG=PRG, MAE=MAE, BIAS=BIAS, noPresencePoints=nrow(inputDF[inputDF$observed ==1,]))
  return(result)
}
