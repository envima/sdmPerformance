#' @title SDM performance evaluation
#'
#' @description
#' Evaluates the performance of species distribution models (SDMs). 
#' The evaluation supports three data types: presence-absence (PA), presence-artificial-absence (PAA), and presence-background (PBG). 
#' Metrics include predictive accuracy (AUC, COR), predictive error measures (MAE, BIAS), across cross-validation replicates.
#'
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
#'   \item{indexPA}{Performance metrics using provided absence data.}
#'   \item{indexPAA}{Performance metrics using artificially generated absence data (AOA-based).}
#'   \item{indexPBG}{Performance metrics using background data.}
#' }
#'
#' @details
#' The index combines accuracy metrics and a spatial correlation stability metric across folds to
#' estimate the reliability and generalizability of the model. COR and AUC measure predictive performance,
#' while layer-wise correlation (Kappa) measures spatial consistency across folds.
#'
#' @examples
#' \dontrun{
#'   result <- evaluateSDMPerformance(
#'               prediction = prediction_raster,
#'               presence = presence_points,
#'               background = TRUE,
#'               environmentalVariables = env_rasters,
#'               replicates = 50)
#'   result$indexPA
#'   result$indexPAA
#'   result$indexPBG
#' }

#'
#' @import terra sf dplyr CAST Metrics prg
#' @export

source("R/run8/functions/generateBackgroundPoints.R")
source("R/run8/functions/generateAAPoints.R")
source("R/run8/functions/calculateMetrics.R")
source("R/run8/functions/indexCalculation.R")

performanceEstimation <- function(
    prediction,
    presence = NA,
    absence = FALSE,
    background = TRUE,
    aa = TRUE,
    environmentalVariables = NA,
    noPointsTesting = NA,
    replicates=100
) {
 
  # -------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------
  
  requireNamespace("sf")
  requireNamespace("terra")
  requireNamespace("CAST")
  requireNamespace("Metrics")
  requireNamespace("dplyr")
  
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
  
  # -------------------------------------------------------------------
  # 1. Presence-Background (PBG)
  # -------------------------------------------------------------------
  if (isTRUE(background)) {
    message(paste("Calculating metrics on presence-background with", replicates, "replicates."))
    indexPBG <- do.call("rbind", lapply(1:replicates, function(i) {
      bg <- generateBackgroundPoints(environmentalVariables, noPointsTesting)
      calculateMetrics(prediction, presence, bg, type = "PBG")
    }))
    indexPBG <- indexPBG %>% dplyr::summarize_all(mean, na.rm = TRUE)
  } else indexPBG <- NA
  
  # -------------------------------------------------------------------
  # 2. Presence-Artificial-Absence (PAA)
  # -------------------------------------------------------------------
  if (is.logical(aa) && isTRUE(aa)) {
    message(paste("Calculating metrics on presence-artificial-absence with", replicates, "replicates."))
    
    extr <- terra::extract(environmentalVariables, presence, ID = FALSE)
    aoa_result <- suppressMessages(CAST::aoa(newdata = environmentalVariables, train = extr, variables = "all", verbose = FALSE))
    aa_mask <- aoa_result$AOA
    aa_mask[aa_mask > 0] <- NA
    
    aa_df <- suppressMessages(as.data.frame(predicts::backgroundSample(aa_mask, n = noPointsTesting * 10, tryf = 30)))
    if(nrow(aa_df) > noPointsTesting) aa_df <- aa_df %>% dplyr::slice_sample(n = noPointsTesting)
    aa <- sf::st_as_sf(aa_df, coords = c("x","y"), crs = terra::crs(aa_mask), remove = FALSE)
    
    indexPAA <- do.call("rbind", lapply(1:replicates, function(i) {
      
      # Sample AA points from precomputed mask
      aa <- generateAAPoints(aa_mask, noPointsTesting)
      
      # Combine presence and artificial absence points
      inputPAA <- na.omit(rbind(
        data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
        data.frame(predicted = terra::extract(prediction, aa)[[2]], observed = 0)
      ))
      
      # Calculate metrics
      indexCalculation(inputPAA, prediction = prediction)
    }))
    indexPAA <- indexPAA %>% dplyr::summarize_all(mean, na.rm = TRUE)
  } else indexPAA <- NA
  
  # -------------------------------------------------------------------
  # 3. Presence-Absence (PA)
  # -------------------------------------------------------------------
  if (!is.logical(absence) || !isFALSE(absence)) {
    indexPA <- calculateMetrics(prediction, presence, absence, type = "PA")
  } else indexPA <- NA
  
  # -------------------------------------------------------------------
  # Combine results
  # -------------------------------------------------------------------
  data <- list(indexPA = indexPA, indexPAA = indexPAA, indexPBG = indexPBG)
  data <- do.call(rbind, Filter(Negate(is.na), data))
  
  gc()
  return(data)
}


