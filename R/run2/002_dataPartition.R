#'@name 001_dataPartition.R
#'@date 25.06.2025
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
library(parallel)
if (Sys.info()[[4]]=="PC19674") setwd("M:/user/bald/SDM/sdmPerformance/")

# 2 - load species ####
#------------------------#

# we picked the VS05 as an example with 120 PA points available
species=list.files("data/PA/",pattern=".gpkg",full.names=F)
#size=c("knndm","random","clusters","block1","block2" ) 

# 2 - split data ####
#--------------------#

if(! file.exists("data/bg.gpkg")){
  bg=as.data.frame(predicts::backgroundSample(terra::rast("data/variables.tif"), n=10000))
  bg=sf::st_as_sf(bg, coords=c("x","y"), crs="epsg:3577", remove=F)
  
  extr=terra::extract(vars,background,ID=F)
  background$Real<- 0
  background$Observed<-0
  background=cbind(background,extr)
  background$random <-  sample(1:6, size=nrow(background), replace = T)
  background$KNNDM <-  sample(1:6, size=nrow(background), replace = T)
  background$clusters <-  sample(1:6, size=nrow(background), replace = T)
  background$block1 <-  sample(1:6, size=nrow(background), replace = T)
  background$block2 <-  sample(1:6, size=nrow(background), replace = T)
  background$x<-NULL
  background$y<-NULL
  
  
  sf::write_sf(bg, "data/bg.gpkg")
} else {
  bg=sf::read_sf("data/bg.gpkg")}

vars=terra::rast("data/variables.tif")


# 3 - create fold in different ways ####
#--------------------------------------#

mclapply(1:length(species), function(i){
  print(i)
  if(!file.exists(paste0("data/run1/virtualSpeciesTrain/",species[i]))){
    vs=sf::read_sf(paste0("data/PA/",species[i]))
    
    # K-nearest neighor distance matching
    KNNDM=CAST::knndm(tpoints=vs,modeldomain = vars,k=6)
    vs$KNNDM <- KNNDM$clusters
    
    # random
    vs$random <- sample(1:6, size=nrow(vs), replace = T)
    
    
    # block cv 1
    block1 = blockCV::cv_spatial(x = vs,
                                 column = "Real",
                                 r = vars, # optionally add a raster layer
                                 k = 6, 
                                 size = 300000, #in m
                                 hexagon = T, 
                                 selection = "random",
                                 offset = c(0, 0),
                                 progress = T, # turn off progress bar for vignette
                                 iteration = 50, 
                                 biomod2 = F,
                                 extend = 5,
                                 plot=F)
    vs$block1 <- block1$folds_ids
    
    
    block2 = blockCV::cv_spatial(x = vs,
                                 column = "Real",
                                 r = vars, # optionally add a raster layer
                                 k = 6, 
                                 size = 100000, #in m
                                 hexagon = F, 
                                 selection = "random",
                                 offset = c(0, 0),
                                 progress = T, # turn off progress bar for vignette
                                 iteration = 50, 
                                 biomod2 = F,
                                 extend = 5,
                                 plot=F)
    vs$block2 <- block2$folds_ids
    
    
    clusters = blockCV::cv_cluster(x=vs,r=vars,k=6)
    vs$clusters <- clusters$folds_ids
    
    # save folding technique
    if(!dir.exists("data/run2/folds")) dir.create("data/run2/folds", recursive=T)
    saveRDS(block1, sprintf("data/run2/folds/%s_block1.RDS",gsub(".gpkg","",species[i])))
    saveRDS(block2, sprintf("data/run2/folds/%s_block2.RDS",gsub(".gpkg","",species[i])))
    saveRDS(clusters, sprintf("data/run2/folds/%s_clusters.RDS",gsub(".gpkg","",species[i])))
    saveRDS(KNNDM, sprintf("data/run2/folds/%s_KNNDM.RDS",gsub(".gpkg","",species[i])))
    
    
    # extract environmental information
    extr=terra::extract(vars,vs,ID=F)
    vs=cbind(vs,extr)
    
    if(!dir.exists("data/run2/virtualSpeciesTrain")) dir.create("data/run2/virtualSpeciesTrain", recursive=T)
    sf::write_sf(vs,paste0("data/run2/virtualSpeciesTrain/",species[i]))
  }
},mc.cores=50)

