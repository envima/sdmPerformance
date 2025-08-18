#' Train a Lasso-Penalized GLM model with spatial cross-validation
#'
#' @description
#' Trains a Lasso-penalized generalized linear model (GLM) using spatially structured 
#' cross-validation and optionally performs raster-based prediction.
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
#' - Uses `caret::train()` with `method = "glmnet"` and `alpha = 1` (Lasso).
#' - Quadratic feature expansion is performed on all predictor variables.
#' - Spatial cross-validation is implemented via `CAST::CreateSpacetimeFolds()`.
#'
#' @examples
#' \dontrun{
#'   lassoModelTraining(df, "presence", c("bio1", "bio12"), "model.rds", "pred.tif", "region", k=5, prediction=TRUE, variables=rasterStack)
#' }
#'
#' @import caret
#' @import CAST
#' @import terra
#' @import climateStability
#' @import assertthat
#' @export



lassoModelTraining <- function(trainingData, # df
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
    savePredictions = "final"#,
    # classProbs = TRUE,
    #summaryFunction = twoClassSummary
  )
  
  # Calculate sample weights to account for class imbalance
  prNum_scv <- as.numeric(table(trainingData[response])["1"])
  bgNum_scv <- as.numeric(table(trainingData[response])["0"])
  
  wt_scv <- ifelse(trainingData[response] == 1, 1, prNum_scv/bgNum_scv) # down-weighting
  
  # clean data
  trainingData<- trainingData%>%dplyr::select(all_of(c(response,predictors)))
  
  
  # generating the quadratic terms for all continuous variables
  # function to creat quadratic terms for lasso and ridge
  quad_obj <- make_quadratic(trainingData, cols = names(variables), verbose = FALSE)
  # now we can predict this quadratic object on the training and testing data
  # this make two columns for each covariates used in the transformation
  training_quad <- predict.make_quadratic(quad_obj, newdata = trainingData)
  testing_quad <- predict.make_quadratic(quad_obj, newdata = raster::stack(variables))
  testing_quad=terra::rast(testing_quad);gc()
  # convert the data.frames to sparse matrices
  
  
  # select all quadratic (and non-quadratic) columns
  new_vars <- setdiff(names(training_quad), response)
  
  
  # --- Lasso grid (alpha = 1, tune over lambda)
  grid <- expand.grid(alpha = 1, lambda = 10^seq(-4, 1, length = 50))
  
  # do classification
  training_quad[[response]] <- as.factor(training_quad[[response]])
  
  # --- Fit Lasso model via caret
  model <- caret::train(
    x = training_quad[new_vars],
    y = training_quad[[response]],
    method = "glmnet",
    trControl = cv_control,
    tuneGrid = grid,
    metric = "Accuracy",
    weights = wt_scv,
    family = "binomial"
  )
  
  
  # If model training fails, return NULL
  if(is.null(model)) return (NULL)
  # Save trained model to disk
  saveRDS(model, outputPathModel)
  
  
  
  # Prediction step (if enabled)
  if(prediction == T){
    pred=terra::predict(testing_quad,model,type="prob")[[2]];rm(testing_quad);gc()
    
    # rescale raster between 0 and 1
    # load original distribution 
    pred=climateStability::rescale0to1(pred)
    names(pred)<-"pred"
    terra::writeRaster(pred, outputPathPrediction,overwrite=T)
  }
}


#' Orthogonal quadratic polynomials for glmnet
#'
#' A function to creat quadratic terms for glmnet functions i.e. lasso and ridge regression.
#' The output is an object of make_quadratic that can be used to predict on rasters and data.frames
#' for creating the quadratic terms.
#'
#' @param df a data.frame, typically the training data.
#' @param cols the name or index of the columns to be transformed. If NULL, all the columns will be transformed.
#' The factor columns won't be transfromed.
#'
#' @author Roozbeh Valavi
#'
#' @return an object of make_quadratic that can be used to predict on rasters and data.frames
#' @export
#'
#' @examples
make_quadratic <- function(df, cols = NULL, verbose = TRUE){
  if(is.null(cols)){
    cols <- colnames(df)
  }
  if(is.numeric(cols)){
    cols <- colnames(df)[cols]
  }
  # remove the factors
  if(any(sapply(df[,cols, drop = FALSE], is.factor))){
    if(verbose){
      message("The factor columns were removed form cols: ", cols[which(sapply(df[,cols, drop = FALSE], is.factor))])
    }
    cols <- cols[-which(sapply(df[,cols, drop = FALSE], is.factor))]
  }
  if(!all(is.element(cols, colnames(df)))){
    stop("The cols should be the same as the column names.")
  }
  xbar <- apply(df[,cols, drop = FALSE], 2, mean)
  x1 <- data.frame(mapply(`-`, df[,cols, drop = FALSE], xbar, SIMPLIFY = FALSE))
  alpha <- colSums(x1 ^ 3) / colSums(x1 ^ 2)
  # specify the output class
  finalList <- list(names = cols, xbars = xbar, alphas = alpha)
  class(finalList) <- c("make_quadratic")
  return(finalList)
}


#' @export
#' @method predict make_quadratic
predict.make_quadratic <- function(object, newdata, cols_from_obj = TRUE, ...){
  if(!methods::is(object, "make_quadratic"))
    stop("object should be a make_quadratic object.")
  # if(!all(object$names %in% names(newdata)))
  #   stop("The newdata does not have the same names as the object.")
  # ncl <- object$names
  # if(!all(names(newdata) %in% object$names))
  #   stop("The newdata does not have the same names as the object.")
  if(cols_from_obj){
    ncl <- object$names
  } else{
    ncl <- names(newdata)
  }
  if(methods::is(newdata, "Raster")){
    for(i in ncl){
      x1 <- newdata[[i]] - object$xbars[i]
      x2 <- (x1 ^ 2) - (object$alphas[i] * x1)
      if(raster::nlayers(newdata) > 1){
        newdata <- newdata[[-which(names(newdata) == i)]]
        newdata <- raster::stack(newdata, x1)
      } else{
        newdata <- x1
      }
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_1")
      newdata <- raster::stack(newdata, x2)
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_2")
    }
  } else if(methods::is(newdata, "data.frame")){
    for(i in ncl){
      x1 <- newdata[,i] - object$xbars[i]
      x2 <- x1 ^ 2 - object$alphas[i] * x1
      newdata[,ncol(newdata) + 1] <- x1
      names(newdata)[ncol(newdata)] <- paste0(i, "_1")
      newdata[,ncol(newdata) + 1] <- x2
      names(newdata)[ncol(newdata)] <- paste0(i, "_2")
      newdata <- newdata[,-which(names(newdata) == i)]
    }
  } else stop("newdata should be a raster or a data.frame.")
  return(newdata)
}

