#' Train a Species Distribution Model Using a Unified Interface
#'
#' @description
#' Trains a species distribution model (SDM) using one of several algorithms, including Random Forest (RF),
#' Maxent, Generalized Additive Models (GAM), Lasso regression, or Boosted Regression Trees (BRT), with spatially structured cross-validation.
#' Model selection is handled via the `modelType` argument, and optional spatial predictions on raster data can be performed.
#'
#' @param trainingData data.frame. Dataset containing species occurrence/background information and predictor variables.
#' @param response character. Name of the binary response variable column (coded as 0 and 1).
#' @param predictors character vector. Names of columns in `trainingData` used as predictor variables.
#' @param outputPathModel character. File path to save the trained model object (e.g., .rds).
#' @param outputPathPrediction character. File path to save the prediction raster (e.g., .tif).
#' @param spacevar character. Column name in `trainingData` defining spatial grouping for cross-validation folds.
#' @param modelType character. One of `"RF"`, `"Maxent"`, `"GAM"`, `"Lasso"`, or `"BRT"` specifying the model algorithm to use.
#' @param k integer. Number of spatial cross-validation folds (default = 5).
#' @param prediction logical. If TRUE, prediction is performed on raster layers provided via `variables`.
#' @param variables SpatRaster or Raster* object. Environmental raster layers used for prediction (required if `prediction = TRUE`).
#' @param xcol character. Name of the column with x-coordinates (used by Maxent only).
#' @param ycol character. Name of the column with y-coordinates (used by Maxent only).
#' @param fc character. Feature class for Maxent tuning (default = `"L"`).
#' @param rm numeric. Regularization multiplier for Maxent tuning (default = 1).
#'
#' @return NULL. The function saves model objects and (optionally) prediction rasters to disk.
#'
#' @details
#' This function acts as a wrapper to streamline training of species distribution models using spatially independent cross-validation.
#' It sources and calls model-specific training functions depending on the `modelType` argument. Each function must be stored in `R/functions/`.
#'
#' The following methods are supported:
#' - `"RF"`: Random Forest
#' - `"Maxent"`: Maximum Entropy Modeling (Phillips et al. 2006)
#' - `"GAM"`: Generalized Additive Models
#' - `"Lasso"`: L1-penalized logistic regression
#' - `"BRT"`: Boosted Regression Trees
#' @examples
#' \dontrun{
#'   trainSpeciesDistributionModel(trainingData = df,
#'                                 response = "presence",
#'                                 predictors = c("bio1", "bio12"),
#'                                 outputPathModel = "model.rds",
#'                                 outputPathPrediction = "prediction.tif",
#'                                 spacevar = "region",
#'                                 modelType = "Maxent",
#'                                 k = 5,
#'                                 prediction = TRUE,
#'                                 variables = rasterStack,
#'                                 xcol = "X", ycol = "Y", fc = "LQH", rm = 1)
#' }
#'
#' @export

trainSpeciesDistributionModel <- function(trainingData,
                                          response,
                                          predictors,
                                          outputPathModel,
                                          outputPathPrediction,
                                          spacevar,
                                          modelType = "RF", # One of: "RF", "Maxent", "GAM", "Lasso", "BRT"
                                          k = 5,
                                          prediction = TRUE,
                                          variables = NULL,
                                          xcol = "X",
                                          ycol = "Y",
                                          fc = "L",
                                          rm = 1) {
  
  # Source individual model functions
  source("R/functions/brtModelTraining.R")
  source("R/functions/gamModelTraining.R")
  source("R/functions/lassoModelTraining.R")
  source("R/functions/maxentModelTraining.R")
  source("R/functions/rfModelTraining.R")
  
  # Dispatch to the appropriate model training function
  if (modelType == "RF") {
    rfModelTraining(trainingData, response, predictors,
                    outputPathModel, outputPathPrediction,
                    spacevar, k, prediction, variables)
    
  } else if (modelType == "Maxent") {
    maxentModelTraining(trainingData, response, predictors,
                        outputPathModel, outputPathPrediction,
                        spacevar, k, prediction, variables,
                        xcol = xcol, ycol = ycol,
                        fc = fc, rm = rm)
    
  } else if (modelType == "GAM") {
    gamModelTraining(trainingData, response, predictors,
                     outputPathModel, outputPathPrediction,
                     spacevar, k, prediction, variables)
    
  } else if (modelType == "Lasso") {
    lassoModelTraining(trainingData, response, predictors,
                       outputPathModel, outputPathPrediction,
                       spacevar, k, prediction, variables)
    
  } else if (modelType == "BRT") {
    brtModelTraining(trainingData, response, predictors,
                     outputPathModel, outputPathPrediction,
                     spacevar, k, prediction, variables)
    
  } else {
    stop("Invalid modelType. Must be one of: 'RF', 'Maxent', 'GAM', 'Lasso', 'BRT'")
  }
}