#'@name 002_modeling.R
#'@date 28.02.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@author 
#'@description train a model for each virtual species
#'@description code taken from the benchmark study of Valavi et al. 2023 DOI: https://doi.org/10.1111/geb.13639
#'@misc OSF repository of the study of Valavi et al. 2023: DOI: 10.17605/OSF.IO/G6DC3
#'@source https://osf.io/puk8v (of the code for the models)

# 1 - install and load packages  ####
#-----------------------------------#



library(tidyverse)
library(sf)
library(parallel)
options(java.parameters = "-Djava.awt.headless=true")
library(rJava)

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=25
}

df=expand.grid(size=as.character(c("KNNDM","random","block1","block2","clusters")) ,
               species=c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
               points=unique(sapply(strsplit(gsub(".gpkg", "", list.files("data/run2/virtualSpeciesTrain", full.names = FALSE)), "_"), `[`, 2)),
               replicates=1:5,
               model=as.character(c("BRT","RF","GAM","Lasso","Maxent")),
               #   model=as.character(c("Lasso")),
               testData = 1:6)

df=df[sample(nrow(df)),]



# 2 - train model ####
#--------------------#

# extract the relevant columns for modelling

source("R/functions/trainSpeciesDistributionModel.R")
source("R/functions/performanceStabilityIndexReplicates.R")

# create a unique name for different runs
nameRun <- paste0("run", 2)


#lapply(1:nrow(df), function(i){
mclapply(1:nrow(df), function(i){
  print(i)
  if(!file.exists(paste0("data/",nameRun,"/maps/",as.character(df$species[i]),"_",df$size[i],"_",df$model[i],"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".tif"))){
    if(!dir.exists(paste0("data/",nameRun,"/models"))) dir.create(paste0("data/",nameRun,"/models"),recursive=T)
    if(!dir.exists(paste0("data/",nameRun,"/maps"))) dir.create(paste0("data/",nameRun,"/maps"),recursive=T)
    
    vars=terra::rast("data/variables.tif")
    
    vs=sf::read_sf(paste0("data/run2/virtualSpeciesTrain/",as.character(df$species[i]),"_",as.character(df$points[i]),"_",df$replicates[i],".gpkg"))
    
    train <- vs %>%
      dplyr::filter(.data[[as.character(df$size[i])]] != df$testData[i]) %>%
      dplyr::filter(Real == 1)
    
    test <- vs %>%
      dplyr::filter(.data[[as.character(df$size[i])]] == df$testData[i])
    
    background=sf::read_sf("data/bg.gpkg")%>%
      dplyr::filter(.data[[as.character(df$size[i])]] != df$testData[i])
    
    
    train<-rbind(train, background)
    
    # get coordinates for maxent modeling
    coords <- as.data.frame(sf::st_coordinates(train))
    train$X <- coords$X
    train$Y<- coords$Y
    train=as.data.frame(train)%>%dplyr::select(-"geom");rm(background,coords)
    
    
    trainSpeciesDistributionModel(trainingData=train, # df
                                  response="Real", # string
                                  predictors=names(vars), # vector with stringds of the predictor columns
                                  outputPathModel=paste0("data/",nameRun,"/models/",as.character(df$species[i]),"_",df$size[i],"_",df$model[i],"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".RDS"), # string
                                  outputPathPrediction=paste0("data/",nameRun,"/maps/",as.character(df$species[i]),"_",df$size[i],"_",df$model[i],"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".tif"), # string of the path to the model output
                                  spacevar=as.character(df$size[i]), # string to the column that holds the assignment of each row to a fold
                                  k=5, # number of folds
                                  prediction=T, # should the predcition also be done?
                                  variables=vars, # if the prediction is done also environmental rasters are needed
                                  modelType = as.character(df$model[i]),
                                  xcol = "X",
                                  ycol = "Y",
                                  fc = "L",
                                  rm = 1
    )
  }
},mc.cores=nCores)
#})

