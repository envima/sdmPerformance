#'@name 005_train_model.R
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
library(terra)
library(sf)
library(parallel)



# read required functions and variables
#source("R/functions/prediction_helper.R")
#source("R/functions/evaluation.R")

# 2 - train model ####
#--------------------#


df=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
               points=c(40,120,200,400),
               replicates=1:10,
               cv="random")
#df$points <- sample(c(40,80,120,160,200,400),size=nrow(df),replace=T)
#species="species034"

#df=df%>%dplyr::filter(df$vs == "VS01")


# get variable names
if(!exists("variables")) variables=terra::rast("data/variables.tif")

mclapply(1:nrow(df), function(s){
#lapply(1:nrow(df), function(s){
  
  
  cat("Calculating model", s, "of", nrow(df), "\n")
  
  bg=sf::read_sf("data/bg.gpkg")
  if(!file.exists(sprintf("data/processed/dataPartition_%s/PA/%s_%s_%s.gpkg",df$cv[s], df$vs[s],df$points[s], df$replicates[s]))) return(NULL)
  sp_presence=read_sf(sprintf("data/processed/dataPartition_%s/PA/%s_%s_%s.gpkg",df$cv[s], df$vs[s],df$points[s], df$replicates[s]))
  
  
  # only presence only data
  sp_presence=sp_presence%>%dplyr::filter(Real==1)
  bg$occ<-0
  sp_presence$occ<-1
  
  # sample background data
  
  if(!file.exists("data/bg_random.gpkg")){
    bg$fold <- sample(1:6, nrow(bg), replace = T)
    environmentalValues=terra::extract(variables, bg, ID=F)
    bg=cbind(bg, environmentalValues)
    sf::write_sf(bg, "data/bg_random.gpkg")
  } else {bg <- sf::read_sf("data/bg_random.gpkg")}
  
  
  
  bg=bg%>%dplyr::select(c(names(variables),"occ", "fold"))
  
  sp_presence=sp_presence%>%dplyr::select(c(names(variables),"occ", "fold"))
  
  # add background data
  pr_bg <- rbind(sp_presence, bg)
  
  pr_bg<-as.data.frame(pr_bg)
  pr_bg <- na.omit(pr_bg)
  # rest rownames for spatial tuning
  #  rownames(pr_bg) <- NULL
  #  rm(bg, sp_presence)
  # extract the relevant columns for modelling
  
  myseed <- 2025
  # train and test sets for spatial cv
  for(z in 1:6){
    
    training_scv <- pr_bg%>%dplyr::filter(fold!=z)%>%dplyr::select(c(names(variables),"occ", "fold"))
    training_scv =na.omit(training_scv);rm(bg,sp_presence)
    
    #*******************************************##
    ## 2.3 - MaxEnt  ####
    #-------------------#
    
    if(!dir.exists("data/output/Maxent/models")) dir.create("data/output/Maxent/models", recursive=T)
    # ptm <- proc.time()
    
    if(!file.exists(paste0(getwd(),"/maxent.jar"))){
      # download maxent.jar file if it does not exist already
      download.file(url="https://osf.io/download/msdy6/", destfile = paste0(getwd(),"/maxent.jar"))
    }
    
    
    for(cvFold in unique(training_scv$fold)){
      
      if(!file.exists(sprintf("data/output/Maxent/models/Maxent_%s_%s_%s_%s_testData%s_cvFold%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z,cvFold))){
      trainingFold <- training_scv[training_scv$fold == cvFold,]
      occurrences <- trainingFold$occ # presence (1s) and background (0s) points
      covariates <- trainingFold[, names(variables)] # predictor covariates
      #set.seed(myseed)
      
      #options(java.parameters = "-Xmx90g")  # Set the Java heap space to 90GB
      
      
      mxnt <- predicts::MaxEnt(x = covariates,
                               p = occurrences,
                              # path = "data/output/Maxent/maxent_files", # path to save maxent files
                               removeDuplicates = FALSE)
      
      saveRDS(mxnt, sprintf("data/output/Maxent/models/Maxent_%s_%s_%s_%s_testData%s_cvFold%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z,cvFold))
      }else mxnt=readRDS(sprintf("data/output/Maxent/models/Maxent_%s_%s_%s_%s_testData%s_cvFold%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z,cvFold))
      # do prediction
      
      if(!file.exists(sprintf("data/output/Maxent/predictions/Maxent_%s_%s_%s_%s_testData%s_cvFold%s.tif",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z,cvFold))){
        
        #options(java.parameters = "-Xmx90g")  # Set the Java heap space to 90GB
        pred <- predicts::predict(object=mxnt, x=variables, args = "outputformat=cloglog")
        names(pred)<-"pred"
        if(!dir.exists(sprintf("data/output/Maxent/predictions/"))) dir.create(sprintf("data/output/Maxent/predictions/"),recursive=T)
        terra::writeRaster(pred, sprintf("data/output/Maxent/predictions/Maxent_%s_%s_%s_%s_testData%s_cvFold%s.tif",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z,cvFold))
      } 
      #t <- proc.time() - ptm
      
      # save processing time
      rm(mxnt, pred, occurrences, covariates);gc()
      
      
      
      
    } # end cv fold loop
    # create final predcition
    if(!file.exists(sprintf("data/output/Maxent/predictions/Maxent_%s_%s_%s_%s_testData%s_mean.tif",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z))){
      p=list.files("data/output/Maxent/predictions/", pattern=sprintf("Maxent_%s_%s_%s_%s_testData%s_cvFold",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z),full.names=T)
      pred=terra::rast(p)
      pred=terra::mean(pred)
      pred=climateStability::rescale0to1(pred)
      terra::writeRaster(pred, sprintf("data/output/Maxent/predictions/Maxent_%s_%s_%s_%s_testData%s_mean.tif",df$vs[s], df$points[s],df$replicates[s],df$cv[2],z));rm(pred,p)
    }
  }
  
  }, mc.cores=45) # end lapply
#})

