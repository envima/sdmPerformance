#' Train a Generalized Additive Model (GAM) with spline terms using spatial cross-validation
#'
#' @description
#' Trains a GAM with spline smoothers for predictors using spatially structured
#' cross-validation implemented via `CAST::CreateSpacetimeFolds()`.
#' Optionally makes raster-based predictions on environmental layers.
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
#' - Fits a GAM with spline terms using `caret::train()` and method `"gamSpline"`.
#' - Uses spatial cross-validation folds defined by `CAST::CreateSpacetimeFolds()`.
#' - Response variable is converted to a factor with levels `"X0"` and `"X1"`.
#' - Predictors are standardized (mean 0, SD 1) based on training data before modeling.
#'
#' @examples
#' \dontrun{
#'   model <- lassoModelTraining(df, "presence", c("bio1", "bio12"),
#'                              "model.rds", "pred.tif", "region", k=5,
#'                              prediction=TRUE, variables=rasterStack)
#' }
#'
#' @import caret
#' @import CAST
#' @import terra
#' @import climateStability
#' @import assertthat
#' @export



gamModelTraining <- function(trainingData, # df
                               response, # string
                               predictors, # vector with strings representing the predictor column names
                               outputPathModel, # string
                               outputPathPrediction, # string of the path to the model output
                               spacevar, # string to the column that holds the assignment of each row to a fold
                               k=5, # number of folds
                               prediction=TRUE, # Should prediction be performed?
                               variables # if the prediction is done also environmental spatRaster are needed
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
  
  
  
  
  # Create spatial cross-validation folds
  cv_folds=CAST::CreateSpacetimeFolds(trainingData, spacevar = spacevar, k=k)
  
  # Define spatial CV settings
  cv_control <- caret::trainControl(
    method = "cv",
    index = cv_folds$index,
    indexOut = cv_folds$indexOut,
    savePredictions = "final",
    classProbs = TRUE
  )
  
  # Calculate sample weights to account for class imbalance
  prNum_scv <- as.numeric(table(trainingData[response])["1"])
  bgNum_scv <- as.numeric(table(trainingData[response])["0"])
  
  wt_scv <- ifelse(trainingData[[response]] == 1, 1, prNum_scv/bgNum_scv) # down-weighting
  
  # clean data
  trainingData<- trainingData%>%dplyr::select(all_of(c(response,predictors)))
  # do classification
  trainingData[[response]] <- as.factor(trainingData[[response]])
  trainingData[[response]] <- factor(trainingData[[response]], levels = c(0, 1), labels = c("X0", "X1"))
  
  testing_variables <- terra::rast()
  # normalising covariates
  for(v in names(variables)){
    #if(v %in% categoricalvars == FALSE){
    # spatial cv covariates
    meaanv <- mean(trainingData[,v])
    sdv <- sd(trainingData[,v])
    trainingData[,v] <- (trainingData[,v] - meaanv) / sdv
    testing_variables[[v]] <- (variables[[v]] - meaanv) / sdv
    #   }
  };gc();rm(meaanv,sdv)
  
  # Fit the GAM with spline terms via caret using method = "gamSpline"
  model <- caret::train(
    x = trainingData[, predictors],
    y = trainingData[[response]],
    method = "gamSpline",
    family = binomial(),
    trControl = cv_control,
    weights = wt_scv,
    metric = "ROC", 
    tuneLength = 5   # number of tuning parameters to try
  )
  
  # If model training fails, return NULL
  if(is.null(model)) return (NULL)
  # Save trained model to disk
  saveRDS(model, outputPathModel)
  
  # Prediction step (if enabled)
  if(prediction == T){
    pred=terra::predict(testing_variables,model,type="prob")[[2]]#;rm(testing_quad);gc()
    
    # rescale raster between 0 and 1
    # load original distribution 
    pred=climateStability::rescale0to1(pred)
    names(pred)<-"pred"
    terra::writeRaster(pred, outputPathPrediction,overwrite=T)
  }
}


