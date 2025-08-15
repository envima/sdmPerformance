#' Train a Maxent Model with Spatial Cross-Validation
#'
#' @description
#' Trains a Maxent species distribution model using presence-background data with spatially structured
#' cross-validation implemented via the `ENMeval` package. Optionally performs raster-based predictions on
#' environmental layers.
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
#' @param xcol character. Column name in `trainingData` indicating the x-coordinate (longitude or easting). Default is "X".
#' @param ycol character. Column name in `trainingData` indicating the y-coordinate (latitude or northing). Default is "Y".
#' @param fc character. Maxent feature class. Defines the type of feature transformation to be used. Default is "L" (linear).
#' @param rm numeric. Regularization multiplier. Controls model complexity (higher values lead to more constrained models). Default is 1.

#'
#' @return NULL. The function saves the trained model and (optionally) prediction raster to disk.
#'
#' @details
#' - Fits a Maxent model using `ENMeval::ENMevaluate()` with user-defined spatial partitions.
#' - Presence and background records must be indicated in the `response` column (1 = presence, 0 = background).
#' - Environmental predictors must be standardized externally if required.
#' - Spatial cross-validation is defined by grouping presence and background records using `spacevar`.
#' - Predictions are scaled to the [0, 1] range using `climateStability::rescale0to1()`.
#' - This function relies on Maxent (Phillips et al. 2006).
#'
#' @references 
#' Muscarella, R., Galante, P. J., Soley‐Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M., & Anderson, R. P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. *Methods in Ecology and Evolution*, 5(11), 1198–1205. https://doi.org/10.1111/2041-210X.12261
#' Kass, J. M., Muscarella, R., Galante, P. J., Boria, R. A., Soley‐Guardia, M., & Anderson, R. P. (2021).ENMeval 2.0: Redesigned for customizable and reproducible modeling of species’ niches and distributions.*Methods in Ecology and Evolution*, 12(9), 1602–1608. https://doi.org/10.1111/2041-210X.13628
#' Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3–4), 231–259. https://doi.org/10.1016/j.ecolmodel.2005.03.026
#'
#'
#' @examples
#' \dontrun{
#'   maxentModelTraining(df, NULL, c("bio1", "bio12"),
#'                       "model.rds", "pred.tif", "region", k=5,
#'                       prediction=TRUE, variables=rasterStack)
#' }
#'
#' @import ENMeval
#' @import terra
#' @import climateStability
#' @import assertthat
#' @import dplyr
#' @export



maxentModelTraining <- function(trainingData, # df
                             response, # string
                             predictors, # vector with strings representing the predictor column names
                             outputPathModel, # string
                             outputPathPrediction, # string of the path to the model output
                             spacevar, # string to the column that holds the assignment of each row to a fold
                             k=5, # number of folds
                             prediction=TRUE, # Should prediction be performed?
                             variables, # if the prediction is done also environmental spatRaster are needed
                             xcol = "X",
                             ycol = "Y",
                             fc = "L",
                             rm = 1
){
  
  # Load required packages
  options(java.parameters = "-Djava.awt.headless=true")
  requireNamespace("ENMeval")
  requireNamespace("terra")
  requireNamespace("climateStability")
  requireNamespace("assertthat")
  requireNamespace("dplyr")  # only if not using dplyr via NAMESPACE
  
  ## Check input validity
  assertthat::assert_that(is.data.frame(trainingData), msg = "trainingData must be a data.frame")
  assertthat::assert_that(response %in% colnames(trainingData), msg = "Response variable not found in trainingData")
  assertthat::assert_that(all(predictors %in% colnames(trainingData)), msg = "Some predictors are not in trainingData")
  assertthat::assert_that(spacevar %in% colnames(trainingData), msg = "spacevar not found in trainingData")
  assertthat::assert_that(xcol %in% colnames(trainingData), msg = "xcol not found in trainingData")
  assertthat::assert_that(ycol %in% colnames(trainingData), msg = "ycol not found in trainingData")
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
  
  # train Maxent model
  # Ensure spatial fold numbers are integers starting from 1
  trainingData[[spacevar]]<- as.integer(factor(trainingData[[spacevar]]))
  # Remove missing values
  trainingData=na.omit(trainingData)
  # Separate presence and background records based on the 'response' column
  occs <- trainingData%>%dplyr::filter(.data[[response]] == 1)
  bg <- trainingData%>%dplyr::filter(.data[[response]] == 0)
  
  # Extract spatial group identifiers
  bg.grp <- bg[[spacevar]]
  occs.grp <- occs[[spacevar]]
  
  # Select only coordinates and predictors
  occs = occs %>% dplyr::select(all_of(c(xcol, ycol, predictors)))
  bg = bg %>% dplyr::select(all_of(c(xcol, ycol, predictors)))
  
  if(!file.exists(paste0(getwd(),"/maxent.jar"))){
    # download maxent.jar file if it does not exist already
    download.file(url="https://osf.io/download/msdy6/", destfile = paste0(getwd(),"/maxent.jar"))
  }
  
  # Train Maxent model using ENMeval with spatial partitions
  model <- ENMeval::ENMevaluate(occs = occs, 
                                bg = bg, 
                                #  envs=variables,
                                algorithm =  "maxent.jar", 
                                partitions = 'user', 
                                user.grp = list(occs.grp=occs.grp,bg.grp=bg.grp),
                                tune.args = list(fc = fc, rm = rm)
  )
  
  # If model training fails, return NULL
  if(is.null(model)) return (NULL)
  # Save trained model to disk
  saveRDS(model, outputPathModel)
  
  # Prediction step (if enabled)
  if(prediction == T){
    pred=terra::predict(model@models[[1]], variables)
    # rescale raster between 0 and 1
    pred=climateStability::rescale0to1(pred)
    names(pred)<-"pred"
    terra::writeRaster(pred, outputPathPrediction,overwrite=T)
  }
}


