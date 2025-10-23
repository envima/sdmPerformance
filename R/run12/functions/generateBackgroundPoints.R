#' @title Generate Background Points for SDM Evaluation
#'
#' @description
#' Generates random background points for presence-background (PBG) evaluation of species distribution models.
#' Points are sampled from the environmental raster layers and returned as an \code{sf} object.
#'
#' @param rasters A \code{terra::SpatRaster} object containing environmental covariates.
#' @param nPoints Integer. Number of background points to generate.
#'
#' @return An \code{sf} object with randomly sampled background points.
#'
#' @export
generateBackgroundPoints <- function(rasters, nPoints) {
  bg_df <- suppressMessages(as.data.frame(predicts::backgroundSample(rasters, n = nPoints * 5)))
  bg_df <- bg_df %>% dplyr::slice_sample(n = nPoints)
  bg <- sf::st_as_sf(bg_df, coords = c("x","y"), crs = terra::crs(rasters), remove = FALSE)
  return(bg)
}

