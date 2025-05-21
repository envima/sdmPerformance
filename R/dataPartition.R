#'@name 003_data_splitting.R
#'@date 28.02.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@description split data into training validation and testing data
#'@misc: zizka paper
#'@misc: enmeval
#'@misc: https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_1.html
#'@misc: https://hannameyer.github.io/CAST/articles/cast02-AOA-tutorial.html


# 1 - set up ####
#---------------#

library(blockCV)
library(predicts)
library(sf)
setwd("M:/user/bald/SDM/sdmPerformance/")

# 2 - split data ####
#--------------------#

if(! file.exists("data/bg.gpkg")){
  bg=as.data.frame(predicts::backgroundSample(terra::rast("data/variables.tif"), n=10000))
  bg=sf::st_as_sf(bg, coords=c("x","y"), crs="epsg:3577", remove=F)
  sf::write_sf(bg, "data/bg.gpkg")
} else {
  bg=sf::read_sf("data/bg.gpkg")}


size=c("knndm","random",50000, 100000,400000, 500000, 600000 ) 


#for ( i in size){
#  if(!(i %in% c("random", "knndmALL","knndmPO","knndmPA","knndmPAA"))){
#    if(!file.exists(sprintf("data/processed/folds/foldPartition_%s.gpkg", i))){
#      #load background data and create one set of folds for all species
#      bg=sf::read_sf("data/raw/background/randomBG_50000.gpkg")
#      #s=species[1]
#      
#      
#      
#      folds = blockCV::cv_spatial(x = bg,
#                                  # column = "Real",
#                                  # r = rasters, # optionally add a raster layer
#                                  k = 6, 
#                                  size = as.numeric(i), #in m
#                                  hexagon = T, 
#                                  selection = "random",
#                                  offset = c(0, 0),
#                                  progress = T, # turn off progress bar for vignette
#                                  iteration = 50, 
#                                  biomod2 = F,
#                                  extend = 5)
#      
#      if(!dir.exists("data/processed/folds")) dir.create("data/processed/folds", recursive=T)
#      saveRDS(folds, sprintf("data/processed/folds/foldPartition_%s.RDS", i))
#      sf::write_sf(folds$blocks,sprintf("data/processed/folds/foldPartition_%s.gpkg", i));rm(bg)
#    } else folds=readRDS(sprintf("data/processed/folds/foldPartition_%s.RDS", i))
#  }

#grid=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
#                 points=c(40,80,120,160,200,400),
#                 replicates=1:10)

grid=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
               points=c(40,120,200,400),
               replicates=2:10)
#species="species034"
i="random"
#grid=grid%>%dplyr::filter(grid$vs == "VS01")

for (s in 1:nrow(grid)){
  print(s) # s="species034"
  # 2.1 - load virtual species data ####
  #------------------------------------#
  
  # if(!file.exists(paste0("data/processed/dataPartition_",i,"/PA/",s,".gpkg"))){
  virtualSpecies=sf::read_sf(paste0("data/PA/",grid$vs[s],"_",grid$points[s],"_",grid$replicates[s],".gpkg"))
  
  if(i =="random"){
    virtualSpecies$fold <- sample(1:6, size=nrow(virtualSpecies), replace = T)
  } 
  
  st_geometry(virtualSpecies) <- "geom"
  
  if(!exists("variables")) variables=terra::rast("data/variables.tif")
  # extract environmental values
  environmentalValues=terra::extract(variables, virtualSpecies, ID=F)
  virtualSpecies=cbind(virtualSpecies, environmentalValues)
  
  # save data 
  if(!dir.exists(sprintf("data/processed/dataPartition_%s/PA", i))) dir.create(sprintf("data/processed/dataPartition_%s/PA", i), recursive=T)
  sf::write_sf(virtualSpecies, paste0("data/processed/dataPartition_",i,"/PA/",grid$vs[s],"_",grid$points[s],"_",grid$replicates[s],".gpkg"))
  
  #clean up
  rm(environmentalValues, virtualSpecies);gc()
}



