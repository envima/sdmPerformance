#' @title Generate Artificial Absence Points for SDM Evaluation
#'
#' @description
#' Generates artificial absence (AA) points for presence-artificial-absence (PAA) evaluation of species distribution models.
#' Uses Area of Applicability (AOA) to identifiy regions that have a high dissimilarity compared ot presence points.
#'
#' @param rasters A \code{terra::SpatRaster} object containing environmental covariates.
#' @param presence An \code{sf} object of known presence locations used for testing.
#' @param nPoints Integer. Number of artificial absence points to generate.
#'
#' @return An \code{sf} object with sampled artificial absence points.
#'
#' @export
generateAAPoints <- function(aa_mask, nPoints) {
  # Sample points from the precomputed AOA mask
  aa_df <- suppressMessages(
    as.data.frame(predicts::backgroundSample(aa_mask, n = nPoints * 10, tryf = 30))
  )
  
  if (nrow(aa_df) > nPoints) aa_df <- aa_df %>% dplyr::slice_sample(n = nPoints)
  
  aa <- sf::st_as_sf(aa_df, coords = c("x", "y"), crs = terra::crs(aa_mask), remove = FALSE)
  
  return(aa)
}
