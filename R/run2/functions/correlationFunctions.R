#install.packages(c("spdep", "sf"))
library(spdep)
library(sf)

corTrainTest<- function(train,test){
  
  coords_train <- st_coordinates(train)
  coords_test <- st_coordinates(test)
  
  # Spatial weights matrix from train to test points (nearest neighbors)
  library(FNN)
  
  # For each train point, find nearest test points within threshold
  if(nrow(test) <= 5) return(NA)
  if(nrow(train) <= 5) return(NA)
  nn <- get.knnx(coords_test, coords_train, k = 5)
  
  # Average test values near each train point
  test_neighbors_mean <- sapply(1:nrow(coords_train), function(i) {
    idx <- nn$nn.index[i, ]
    mean(test$bio_1[idx], na.rm = TRUE)
  })
  
  # Correlate train values with mean test neighbors values
  correlation <- cor(train$bio_1, test_neighbors_mean, use = "complete.obs")
  return(correlation)
  
}

##################################

moran <- function(test){
  coords <- st_coordinates(test)
  
  coords <- st_coordinates(test)
  dist_matrix <- as.matrix(dist(coords))
  
  summary(dist_matrix[lower.tri(dist_matrix)])[5]
  
  d_threshold <- summary(dist_matrix[lower.tri(dist_matrix)])[5][[1]]  # adjust based on your spatial resolution
  
  # Define neighbors within threshold distance
  nb_test <- dnearneigh(coords, 0, d_threshold)
  
  # Create spatial weights
  lw_test <- tryCatch(
    nb2listw(nb_test, style = "W"),
    error = function(e) {
      message("Spatial weights creation failed: ", e$message)
      return(NULL)
    }
  )
  if(is.null(lw_test)) return(NA)
  
  moran_result <- tryCatch(
     moran.test(test$bio_1, lw_test),
    error = function(e) {
      message("Moran calculation failed: ", e$message)
      return(NULL)
    }
  )
  if(is.null( moran_result)) return(NA)
  
  
  return(moran_result$estimate[1][1])
}
