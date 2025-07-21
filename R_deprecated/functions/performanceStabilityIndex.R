#' @title Performance Stability Index Calculation
#'
#' @description
#' Computes a composite metric quantifying the performance and stability of spatial prediction models
#' using a combination of predictive accuracy (AUC, COR) and inter-layer correlation (spatial stability)
#' across cross-validation folds. The index is derived using the `Lithics3D:::getTArea()` method on a 
#' standardized vector representation of metrics.
#'
#' @param x A `terra::SpatRaster` object (multi-layer raster stack) representing prediction outputs across CV folds.
#' @param absence Optional. A `sf` object representing known absence locations. If not provided, `background` or `aa` must be used.
#' @param presence Required. A `sf` object representing known presence locations.
#' @param background Optional. Either a `sf` object with background points, or logical `TRUE` to auto-generate them. Logical `FALSE` to not calculate them
#' @param aa Optional. A `sf` object representing artificial absence points. If not provided, they are derived using AOA (Area of Applicability).
#' @param environmentalVariables Optional unless `background = TRUE` or `aa` is missing. A `terra::SpatRaster` with environmental covariates.
#' @param noPointsTesting Integer. Number of background or artificial absence points to generate.
#' @param trainingData Presence only data on which the model was trained. A `sf` object.
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
#' @import terra sf dplyr CAST Metrics
#' @export

performanceStabilityIndex <- function(
    x = NA,
    presence = NA,
    absence = FALSE,
    background = TRUE,
    aa = TRUE,
    environmentalVariables = NA,
    trainingData = NA,
    noPointsTesting = NA
) {
  requireNamespace("sf")
  requireNamespace("terra")
  requireNamespace("CAST")
  requireNamespace("Metrics")
  requireNamespace("dplyr")
  
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra::SpatRaster object.")
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
  
  if ((isTRUE(background) ||  isTRUE(aa)) && !inherits(x, "SpatRaster")) {
    stop("Environmental variables must be provided to generate background or artificial absence data.")
  }
  
  if (is.logical(aa) && isTRUE(aa) && is.na(trainingData)[1]) {
    stop("Training data must be provided to calculate artificial absence (AA) points.")
  }
  
  if (is.na(noPointsTesting)) {
    noPointsTesting <- nrow(presence)
  }
  
  pred <- terra::mean(x)
  stability <- stabilityRasters(x)
  
  if (isTRUE(background)) {
    bg_df <- as.data.frame(predicts::backgroundSample(environmentalVariables, n = noPointsTesting*5))
    bg_df <- bg_df%>%dplyr::slice_sample( n=noPointsTesting)
    bg <- sf::st_as_sf(bg_df, coords = c("x", "y"), crs = terra::crs(environmentalVariables), remove = FALSE)
    
    inputPBG <- na.omit(rbind(
      data.frame(predicted = terra::extract(pred, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(pred, dplyr::sample_n(bg, nrow(presence)))[[2]], observed = 0)
    ))
    indexPBG <- indexCalculation(inputPBG, stability)
  } else indexPBG <- NA
  
  if (is.logical(aa) && isTRUE(aa)) {
    extr <- terra::extract(environmentalVariables, trainingData, ID = FALSE)
    aoa_result <- suppressMessages(CAST::aoa(newdata = environmentalVariables, train = extr, variables = "all",verbose = FALSE))
    aa_mask <- aoa_result$AOA
    aa_mask[aa_mask > 0] <- NA
    aa_df <- as.data.frame(predicts::backgroundSample(aa_mask, n = noPointsTesting*5, tryf = 5))
    aa_df <- aa_df%>%dplyr::slice_sample(n=noPointsTesting)
    aa <- sf::st_as_sf(aa_df, coords = c("x", "y"), crs = terra::crs(aa_mask), remove = FALSE)
    
    inputPAA <- na.omit(rbind(
      data.frame(predicted = terra::extract(pred, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(pred, dplyr::sample_n(aa, nrow(presence)))[[2]], observed = 0)
    ))
    indexPAA <- indexCalculation(inputPAA, stability)
  } else indexPAA <- NA
  
  if (!is.logical(absence) || !isFALSE(absence)) {
    inputPA <- na.omit(rbind(
      data.frame(predicted = terra::extract(pred, presence)[[2]], observed = 1),
      data.frame(predicted = terra::extract(pred, absence)[[2]], observed = 0)
    ))
    indexPA <- indexCalculation(inputPA, stability)
  } else indexPA <- NA
  
  list(
    indexPA = indexPA,
    indexPAA = indexPAA,
    indexPBG = indexPBG
  )
}

stabilityRasters <- function(r) {
  corr_matrix <- terra::layerCor(r, fun = "cor")$correlation
  corr_matrix[!lower.tri(corr_matrix, diag = FALSE)] <- NA
  valid_corrs <- na.omit(as.vector(corr_matrix))
  mean(valid_corrs)
}

indexCalculation <- function(inputDF, stability) {
  COR <- if (length(unique(inputDF$predicted)) > 1) cor(inputDF$observed, inputDF$predicted) else NA
  AUC <- Metrics::auc(inputDF$observed, inputDF$predicted)
  values <- c(COR, 0, 0, 0, stability, 0, 0, 0, AUC)
  metric <- Lithics3D:::getTArea(values)
  metric <- metric / 0.8660254
  
  metric=data.frame(metric=metric, AUC=AUC, COR=COR,stability=stability)
  return(metric)
}
