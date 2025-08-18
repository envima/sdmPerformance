#' Train a Random Forest (RF) model with spatial cross-validation
#'
#' @description
#' Trains a Random Forest (RF) model using spatially structured 
#' cross-validation folds and optionally performs raster-based prediction.
#' The function is designed for binary response variables and accounts for class imbalance.
#'
#' @param trainingData data.frame. Tabular dataset containing the response and predictor variables.
#' @param response character. Name of the binary response variable (coded as 0 and 1).
#' @param predictors character vector. Names of the predictor variable columns.
#' @param outputPathModel character. File path for saving the trained model (RDS format).
#' @param outputPathPrediction character. File path for saving the raster predictions (GeoTIFF or other format supported by terra).
#' @param spacevar character. Column name in `trainingData` defining spatial grouping for cross-validation folds.
#' @param k integer. Number of cross-validation folds (default is 5).
#' @param prediction logical. If TRUE, model predictions will be made on the provided raster data (`variables`).
#' @param variables SpatRaster or Raster* object. Environmental layers used for prediction (required if `prediction = TRUE`).
#'
#' @return NULL. The function saves the trained model and (optionally) prediction raster to disk.
#'
#' @details
#' - Uses the `caret` package for model training.
#' - Spatial folds are generated using `CAST::CreateSpacetimeFolds`.
#' - A fixed grid search is used for hyperparameter tuning.
#' - Random forest is downwsampled, same amount of background points as presence points are used.
#'
#' @examples
#' \dontrun{
#'   rfModelTraining(df, "presence", c("bio1", "bio12"), "model.rds", "pred.tif", "region", k=5, prediction=TRUE, variables=rasterStack)
#' }
#'
#' @import caret
#' @import CAST
#' @import terra
#' @import climateStability
#' @import assertthat
#' @export



rfModelTraining <- function(trainingData, # df
                             response, # string
                             predictors, # vector with stringds of the predictor columns
                             outputPathModel, # string
                             outputPathPrediction, # string of the path to the model output
                             spacevar, # string to the column that holds the assignment of each row to a fold
                             k=5, # number of folds
                             prediction=T, # should the predcition also be done?
                             variables # if the prediction is done also environmental rasters are needed
){
  
  # Load required packages
  requireNamespace("caret")
  requireNamespace("CAST")
  requireNamespace("terra")
  requireNamespace("climateStability")
  requireNamespace("assertthat")
  
  ## Check input validity
  assertthat::assert_that(is.data.frame(trainingData), msg = "trainingData must be a data.frame")
  assertthat::assert_that(response %in% colnames(trainingData), msg = "Response variable not found in trainingData")
  assertthat::assert_that(all(predictors %in% colnames(trainingData)), msg = "Some predictors are not in trainingData")
  assertthat::assert_that(spacevar %in% colnames(trainingData), msg = "spacevar not found in trainingData")
  assertthat::assert_that(is.character(response), is.character(predictors), is.character(spacevar))
  assertthat::assert_that(is.character(outputPathModel), is.character(outputPathPrediction))
  assertthat::assert_that(is.numeric(k) && k > 1)
  assertthat::assert_that(is.logical(prediction))
  if (prediction) {
    assertthat::assert_that(inherits(variables, "SpatRaster") || inherits(variables, "Raster"), 
                            msg = "variables must be a SpatRaster or Raster object")
  }
  
  # Ensure binary response is coded as 0/1
  response_vals <- unique(trainingData[[response]])
  assertthat::assert_that(all(response_vals %in% c(0, 1)), msg = "Response must be binary with values 0 and 1")
  
  s=min(table(trainingData[[response]]))
  trainingData<-  trainingData %>%
    group_by(.data[[response]]) %>%
    sample_n(s) %>%
    ungroup()
  
  trainingData[[response]] <- as.factor(trainingData[[response]])
  
  # Create spatial cross-validation folds
  cv_folds=CAST::CreateSpacetimeFolds(trainingData, spacevar = spacevar, k=k)
  
  # Define spatial CV settings
  cv_control <- caret::trainControl(
    method = "cv",
    index = cv_folds$index,
    indexOut = cv_folds$indexOut,
    savePredictions = "final"
  )
  
  # Tune grid for Random Forest
  grid <- expand.grid(
    mtry = 1:length(predictors)) 
  
  
  # Train the Random Forest model using caret
  RF <- caret::train(
    x = trainingData[predictors],
    y = trainingData[[response]],
    method = "rf",
    trControl = cv_control,
    tuneGrid = grid,
    ntree = 1000,
    importance = TRUE
  )
  
  
  # If model training fails, return NULL
  if(is.null(RF)) return (NULL)
  # Save trained model to disk
  saveRDS(RF, outputPathModel)
  
  # Prediction step (if enabled)
  if(prediction == T){
    pred=terra::predict(variables, RF, type = "prob",na.rm=T)[[2]]
    names(pred)<-"pred"
    # rescale raster between 0 and 1
    pred=climateStability::rescale0to1(pred)
    terra::writeRaster(pred, outputPathPrediction,overwrite=T)
  }
  
}

