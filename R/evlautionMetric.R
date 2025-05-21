#R function to evaluate SDMS
library(dplyr)
library(Lithics3D)
# dependencies: terra, sf, Lithics3D, predics, CAST

#x=replcate rasters as spatRasters
#Absence
#Presence 
#samppleBackground
#background
#PAA
#metrics= AUC, 
setwd("M:/user/bald/SDM/sdmPerformance/")
p=list.files("data/output/Maxent/predictions/", pattern="Maxent_VS01_40_1_random_testData4_cvFold", full.names=T)
pred=terra::rast(p)
pa=sf::read_sf("data/processed/dataPartition_random/PA/VS01_40_1.gpkg")%>%dplyr::filter(fold == 1)
presence=pa%>%dplyr::filter(Real==1)
absence=pa%>%dplyr::filter(Real==0)
trainingData=sf::read_sf("data/processed/dataPartition_random/PA/VS01_40_1.gpkg")%>%dplyr::filter(fold != 1)%>%dplyr::filter(Real==1)

rm(pa,pred,p)

performanceStabilityIndex<- function(
    x, # raster stack with cv folds
    absence = NA,# or sf object
    presence = NA, 
    background= NA ,# (or dataframe with points) or TRUE 
    aa = NA,
    environmentalVariables=NA, # terra::rast() object
    noBackgroundPoints=10000
){
  
  #if (is.na(absence) && background == F %% artificalAbsence == F){
  #  return("Either absence, background or aa points are needed. Do either of the below: turn on of the paramters vackground or artificial absence to true or provdie a dataset for the parameters absence, background or aa in sf format.")
  #}
  
  #if((background == T | aa == T) && is.na(environmentalVariables)){
  #  return("to calculate backgorund points and / or artificalAbsenc epoints a raster of the environmental vairables needs to be given")
  #}
  
  
  # 1 - create background and artificial Absence datasets
  
  bg=as.data.frame(predicts::backgroundSample(environmentalVariables, n=10000))
  bg=sf::st_as_sf(bg, coords=c("x","y"), crs=terra::crs(environmentalVariables), remove=F)
  
  
  extr=terra::extract(environmentalVariables, trainingData,ID=F)
  
  
  AOA=CAST::aoa(newdata=environmentalVariables,
                train = extr,
                variables = "all")
  aa=AOA$AOA
  aa[aa > 0 ]<-NA
  aa=as.data.frame(predicts::backgroundSample(aa, n=10000, tryf=5))
  aa=sf::st_as_sf(aa, coords=c("x","y"), crs=terra::crs(AOA$AOA), remove=F)
  
  
  
  
  # 2 - calculate metrics ####
  pred=terra::mean(x)
  
  inputPA = na.omit(rbind(data.frame(predicted=terra::extract(pred,presence)[[2]],
                                     observed=1),
                          data.frame(predicted=terra::extract(pred,absence)[[2]],
                                     observed=0)))
  
  indexPA <- indexCalculation(inputPA,x)
  
  
  
  
  inputPAA = na.omit(rbind(data.frame(predicted=terra::extract(pred,presence)[[2]],
                                      observed=1),
                           data.frame(predicted=terra::extract(pred,aa%>%dplyr::sample_n(nrow(presence)))[[2]],
                                      observed=0)))
  indexPAA=indexCalculation(inputPAA,x)
  
  inputPBG = na.omit(rbind(data.frame(predicted=terra::extract(pred,presence)[[2]],
                                      observed=1),
                           data.frame(predicted=terra::extract(pred,bg%>%dplyr::sample_n(nrow(presence)))[[2]],
                                      observed=0)))
  indexPBG<-indexCalculation(inputPBG,x)
  
  # create object to return
  metric@indexPBG <- indexPBG
  
}


stabililtyRasters<- function(r){
  stability=terra::layerCor(r, fun="cor")$correlation
  stability[!lower.tri(stability, diag = FALSE)] <- NA
  stability <- na.omit(as.vector(stability))
  stability<- mean(stability)
  return(stability)
}

indexCalculation<- function(inputDF, r){
  unbalanced=data.frame(COR=cor( inputDF$observed, inputDF$predicted),
                        AUC=Metrics::auc(inputDF$observed,inputDF$predicted),
                        Kappa=stabililtyRasters(r))
  
  # calculate metric
  values=c(unbalanced$COR,0,0,0,unbalanced$Kappa,0,0,0,unbalanced$AUC)
  metric=Lithics3D:::getTArea(values)
  metric= metric / 0.8660254
  return(metric)
}

#grid=expand.grid(stability=seq(0,1,0.1),
#            AUC=seq(0,1,0.1),
#            COR=seq(0,1,0.1),
#            index=NA)#

#for (i in 1:nrow(grid)){
#  # calculate metric
#  values=c(grid$COR[i],0,0,0,grid$stability[i],0,0,0,grid$AUC[i])
#  metric=Lithics3D:::getTArea(values)
#  grid$index[i]<-metric
#}



