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

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

df=expand.grid(size=as.character(c("KNNDM","random","block1","block2","clusters")) ,
               species=c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
               points=unique(sapply(strsplit(gsub(".gpkg", "", list.files("data/run2/virtualSpeciesTrain", full.names = FALSE)), "_"), `[`, 2)),
               replicates=1:5,
               # model=as.character(c("RF")),
               model=as.character(c("Lasso", "RF","Maxent", "BRT", "GAM")),
               testData = 1:6)




# create a unique name for different runs
nameRun <- paste0("run", 2)


# 4 - evaluate models ####
#------------------------#

vars_path=normalizePath("data/variables.tif")


source("R/functions/correlationFunctions.R")
source("R/functions/performanceStabilityIndexReplicates.R")

# Pre-filter already processed files
#df <- df[!file.exists(paste0("data/",nameRun,"/results/",
#                             df$species, "_", df$size, "_", df$model,
#                             "_testData", df$testData,
#                             "_points", df$points,
#                             "_replicates", df$replicates, ".RDS")), ]



mclapply(1:nrow(df), function(i){
  
  print(i)
  if(!file.exists(paste0("data/",nameRun,"/results/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".RDS"))){
    
    vars=terra::rast(vars_path)
    vs=sf::read_sf(paste0("data/run2/virtualSpeciesTrain/",as.character(df$species[i]),"_",as.character(df$points[i]),"_",df$replicates[i],".gpkg"))
    
    train <- vs %>%
      dplyr::filter(.data[[as.character(df$size[i])]] != df$testData[i]) %>%
      dplyr::filter(Real == 1)
    
    test <- vs %>%
      dplyr::filter(.data[[as.character(df$size[i])]] == df$testData[i])
    
    if(nrow(test%>%dplyr::filter(Real==1))<1) return(NULL)
    
    COR<- corTrainTest(train,test)
    moranTrain<- moran(train)
    moranTest<- moran(test)
    
    realDistirbution=terra::rast(paste0("data/virtualSpecies/",as.character(df$species[i]),".tif"))
    pred=terra::rast(paste0("data/",nameRun,"/maps/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".tif"))
    if(isTRUE(terra::global(pred, fun = function(x) all(is.na(x)))[[1]])) return(NULL)
    pred=terra::mask(pred,realDistirbution)
    x=list.files(paste0("data/",nameRun,"/maps/"),pattern=paste0(as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData"),full.names = T)
    pattern <- paste0("points",df$points[i],"_replicates",df$replicates[i],".tif")
    x <- lapply(x, function(vec) vec[grepl(pattern, vec)])
    # Remove empty elements (length 0)
    x <- Filter(function(a) length(a) > 0, x)
    x=terra::rast(unlist(x))
    
    result=performanceStabilityIndex( x = x,
                                      prediction=pred,
                                      presence = test%>%dplyr::filter(Real==1),
                                      absence = test%>%dplyr::filter(Real==0),
                                      background = TRUE,
                                      aa = TRUE,
                                      environmentalVariables = vars,
                                      trainingData = train%>%dplyr::filter(Real==1),
                                      noPointsTesting = NA,
                                      replicates=100)
    
    result$trueCor <- terra::layerCor(terra::rast(list(pred,realDistirbution)),fun="cor")$correlation[[1,2]]
    result$model <- as.character(df$model[i])
    result$size <- as.character(df$size[i])
    result$testData <- df$testData[i]
    result$method <- gsub("index","",rownames(result))
    result$corTrainTest<-COR
    result$moranTrain<-moranTrain
    result$moranTest <- moranTest
    result$replicate <- as.character(df$replicates[i])
    result$points <- as.character(df$points[i])
    result$species <- as.character(df$species[i])
    
    rm(pred,realDistirbution,x,train,test,moranTest,moranTrain,pattern,vs,COR);gc()
    if(!dir.exists(paste0("data/",nameRun,"/results"))) dir.create(paste0("data/",nameRun,"/results"))
    saveRDS(result, paste0("data/",nameRun,"/results/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData",df$testData[i],"_points",as.character(df$points[i]),"_replicates",df$replicates[i],".RDS"))
  }
},mc.cores=nCores)


gc()


if(!file.exists("data/run2/results.RDS")){
  data=list.files("data/run2/results",full.names = T)
  data=mclapply(data, function(x){
    df=readRDS(x)
    df$species <- strsplit(strsplit(x, split="/")[[1]][4],split="_")[[1]][1]
    return(df)
  },mc.cores=nCores)
  
  #data2=do.call(rbind,data)
  data2=dplyr::bind_rows(data)
  saveRDS(data2, "data/run2/results.RDS")
} else data=readRDS("data/run2/results.RDS")


