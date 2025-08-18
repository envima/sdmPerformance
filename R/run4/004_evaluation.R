#'@name 004_evaluation.R
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
library(future)
library(future.apply)
library(spdep)
library(FNN)
library(parallel)
library(RandomFields)
library(RandomFieldsUtils)
library(NLMR)
library(raster)

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

df=expand.grid(size=as.character(c("KNNDM","random","block1","block2","clusters")) ,
               species=c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
               points=unique(sapply(strsplit(gsub(".gpkg", "", list.files("data/run2/virtualSpeciesTrain", full.names = FALSE,pattern=".gpkg")), "_"), `[`, 2)),
               replicates=1:5,
               # model=as.character(c("RF")),
               model=as.character(c("Lasso")),
               testData = 1:6)

set.seed(123)  # optional, for reproducibility
df <- df[sample(nrow(df)), ]


# create a unique name for different runs
nameRun <- paste0("run", 4)


# 4 - evaluate models ####
#------------------------#

# specify RandomFields options ----
RandomFields::RFoptions(cPrintlevel = 0)
RandomFields::RFoptions(spConform = FALSE)
RandomFields::RFoptions(install="no")
# set RF seed ----
RandomFields::RFoptions(seed = NULL)

vars_path=normalizePath("data/variables.tif")

source(paste0("R/",nameRun,"/functions/nlm_gaussianfield.R"))
source(paste0("R/",nameRun,"/functions/RFsimulate_custom.R"))
source(paste0("R/",nameRun,"/functions/correlationFunctions.R"))
source(paste0("R/",nameRun,"/functions/performanceStabilityIndexReplicates.R"))

# Pre-filter already processed files
#df <- df[!file.exists(paste0("data/",nameRun,"/results/",
#                             df$species, "_", df$size, "_", df$model,
#                             "_testData", df$testData,
#                             "_points", df$points,
#                             "_replicates", df$replicates, ".RDS")), ]



mclapply(1:100, function(i){
  
  print(i)
  if(!file.exists(paste0("data/",nameRun,"/results/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".RDS"))){
    
    vars=terra::rast(vars_path)
    vs=sf::read_sf(paste0("data/run2/virtualSpeciesTrain/",as.character(df$species[i]),"_",as.character(df$points[i]),"_",df$replicates[i],".gpkg"))
    
    test <- vs %>%
      dplyr::filter(.data[[as.character(df$size[i])]] == df$testData[i])
    
    if(nrow(test%>%dplyr::filter(Real==1))<1) return(NULL)
    
    
    realDistribution=terra::rast(paste0("data/virtualSpecies/",as.character(df$species[i]),".tif"))
    
    
    #  random values between 0 and 1
    autocorrRange=sample(c(20,50,100,500,800),size=1)
    randomField=nlm_gaussianfield(nrow=1777, ncol=2247, resolution=813.3488,
                                        autocorr_range = autocorrRange)
    
  
    randomField=terra::rast(randomField)
    terra::ext(randomField) <- terra::ext(realDistribution)
    terra::crs(randomField)<- terra::crs(realDistribution)
    randomField=terra::mask(randomField, realDistribution)
    #terra::plot(randomField)
    
    # Weighted
    randomEffects <- sample(seq(0,1,0.05), size=1)
    pred <- (randomEffects * realDistribution) + ((1-randomEffects) * randomField)
    pred=climateStability::rescale0to1(pred)
    
    if(isTRUE(terra::global(pred, fun = function(x) all(is.na(x)))[[1]])) return(NULL)
    
    result=performanceStabilityIndex( x = NA,
                                      prediction=pred,
                                      presence = test%>%dplyr::filter(Real==1),
                                      absence = test%>%dplyr::filter(Real==0),
                                      background = TRUE,
                                      aa = TRUE,
                                      environmentalVariables = vars,
                                      noPointsTesting = NA,
                                      replicates=100)
    
    result$trueCor <- terra::layerCor(terra::rast(list(pred,realDistribution)),fun="cor")$correlation[[1,2]]
    result$model <- as.character(df$model[i])
    result$size <- as.character(df$size[i])
    result$testData <- df$testData[i]
    result$method <- gsub("index","",rownames(result))
    result$corTrainTest<-NA
    result$moranTrain<-NA
    result$moranTest <- NA
    result$replicate <- as.character(df$replicates[i])
    result$points <- as.character(df$points[i])
    result$species <- as.character(df$species[i])
    result$randomEffects <- randomEffects
    result$autocorrRange <- autocorrRange
    
    rm(pred,realDistribution,test,vs, autocorrRange, randomEffects);gc()
    if(!dir.exists(paste0("data/",nameRun,"/results"))) dir.create(paste0("data/",nameRun,"/results"), recursive=T)
    saveRDS(result, paste0("data/",nameRun,"/results/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".RDS"))
  }
},mc.cores=nCores)


gc()


if(!file.exists(paste0("data/",nameRun,"/results.RDS"))){
  data=list.files(paste0("data/",nameRun,"/results"),full.names = T)
  data=mclapply(data, function(x){
    df=readRDS(x)
    df$species <- strsplit(strsplit(x, split="/")[[1]][4],split="_")[[1]][1]
    return(df)
  },mc.cores=nCores)
  
  #data2=do.call(rbind,data)
  data2=dplyr::bind_rows(data)
  saveRDS(data2, paste0("data/",nameRun,"/results.RDS"))
} else data=readRDS(paste0("data/",nameRun,"/results.RDS"))


