source("R/functions/performanceStabilityIndex2.R")
# Libraries
library(tidyverse)
library(parallel)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(tidyr)
library(ggplot2)
library(caret)
library(ggpmisc)

if (Sys.info()[[4]]=="PC19674") setwd("M:/user/bald/SDM/sdmPerformance/")

df=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
               points=c(40,120,200,400),
               replicates=1:10,
               cv="random",
               testData=1:6,
               method=NA,
               index=NA,
               AUC=NA,
               COR=NA,
               stability=NA)


#df=df%>%dplyr::filter(df$vs == "VS01")

if(!dir.exists("data/resultsAPAA")) dir.create("data/resultsAPAA")

mclapply(1:nrow(df), function(s){
  #lapply(1:nrow(df), function(s){
  print(s)
  if(!file.exists(sprintf("data/resultsAPAA/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))){
    df2=df[s,]
    
    # s=1 
    p=list.files("data/output/Maxent/predictions/", pattern=sprintf("Maxent_%s_%s_%s_%s_testData%s_cvFold",df$vs[s], df$points[s],df$replicates[s],df$cv[s],df$testData[s]),full.names=T)
    p=terra::rast(p)
    vs=terra::rast(paste0("data/virtualSpecies/",df$vs[s],".tif"))
    df2$trueCOR<- terra::layerCor(c(vs, terra::mean(p)),fun="cor")$correlation[1,2]
    index=performanceStabilityIndex(x=p,
                                    absence = sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==0)%>%dplyr::filter(fold == df$testData[s]),
                                    presence=sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==1)%>%dplyr::filter(fold == df$testData[s]),
                                    background=TRUE,
                                    trainingData=sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==1)%>%dplyr::filter(fold != df$testData[s]),
                                    aa=TRUE,
                                    environmentalVariables=terra::rast("data/variables.tif"),
                                    noPointsTesting=NA,
                                    artificialPresence=T)
    dfPAA=df2
    dfPAA$index<- index$indexPAA$metric
    dfPAA$AUC<- index$indexPAA$AUC
    dfPAA$COR<- index$indexPAA$COR
    dfPAA$PRG<- index$indexPAA$PRG
    dfPAA$MAE<- index$indexPAA$MAE
    dfPAA$BIAS<- index$indexPAA$BIAS
    dfPAA$X<- index$indexPAA$X
    dfPAA$indexA<- index$indexPAA$indexA
    dfPAA$stability<- index$indexPAA$stability
    dfPAA$method <- "PAA"
    dfPAA$noPresencePoints <- index$indexPAA$noPresencePoints
    
    if (!is.na(index$indexAPAA)[1]){
      dfAPAA=df2
      dfAPAA$index<- index$indexAPAA$metric
      dfAPAA$AUC<- index$indexAPAA$AUC
      dfAPAA$COR<- index$indexAPAA$COR
      dfAPAA$PRG<- index$indexAPAA$PRG
      dfAPAA$MAE<- index$indexAPAA$MAE
      dfAPAA$BIAS<- index$indexAPAA$BIAS
      dfAPAA$X<- index$indexAPAA$X
      dfAPAA$indexA<- index$indexAPAA$indexA
      dfAPAA$stability<- index$indexAPAA$stability
      dfAPAA$method <- "APAA"
      dfAPAA$noPresencePoints <- index$indexAPAA$noPresencePoints
    }
    
    if (!is.na(index$indexPA)[1]){
      dfPA=df2
      dfPA$index<- index$indexPA$metric
      dfPA$AUC<- index$indexPA$AUC
      dfPA$COR<- index$indexPA$COR
      dfPA$stability<- index$indexPA$stability
      dfPA$PRG<- index$indexPA$PRG
      dfPA$MAE<- index$indexPA$MAE
      dfPA$BIAS<- index$indexPA$BIAS
      dfPA$method <- "PA"
      dfPA$X<- index$indexPA$X
      dfPA$indexA<- index$indexPA$indexA
      dfPA$noPresencePoints <- index$indexPA$noPresencePoints
    }
    dfPBG=df2
    dfPBG$index<- index$indexPBG$metric
    dfPBG$X<- index$indexPBG$X
    dfPBG$indexA<- index$indexPBG$indexA
    dfPBG$AUC<- index$indexPBG$AUC
    dfPBG$COR<- index$indexPBG$COR
    dfPBG$stability<- index$indexPBG$stability
    dfPBG$PRG<- index$indexPBG$PRG
    dfPBG$MAE<- index$indexPBG$MAE
    dfPBG$BIAS<- index$indexPBG$BIAS
    dfPBG$method <- "PBG"
    dfPBG$noPresencePoints <- index$indexPBG$noPresencePoints
    
    if (is.na(index$indexPA)[1]){
      df2=rbind( dfPAA, dfPBG,dfAPAA)
    } else {df2=rbind(dfPA, dfPAA, dfPBG,dfAPAA)}
    
    rm(index, p,vs, dfPA, dfPAA, dfPBG,dfAPAA);gc()
    saveRDS(df2, sprintf("data/resultsAPAA/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))
  }
  #})
}, mc.cores=45)


l=list.files("data/resultsAPAA/", full.names=T)
l=do.call(rbind, lapply(l, function(x){
  data=readRDS(x)
  return(data)}))


l$index <- round(l$index,4)
l$diffIndex <- round(abs(l$index - l$trueCOR),4)
l$diffIndex2 <- round(abs(l$index2 - l$trueCOR),4)
l$trueCOR <- abs(l$trueCOR)
#l[l < 0] <- 0
#l=l%>%dplyr::filter(method=="PAA")

# Mean absoulte error of outliers
#observed=l$trueCOR
#predicted=l$COR
#residuals <- observed - predicted
#threshold <- 2 * sd(residuals)
#outliers <- abs(residuals) > threshold
#mae_outliers <- mean(abs(residuals[outliers]))


indexCalculation3 <- function(COR,stability,AUC,MAE,PRG,BIAS, index) {
  
  if(AUC < 0.5) AUC<-0.5
  aucRescale=data.frame(AUCnew=c(rep(0,21),seq(0,1,0.05)),
                        AUC=c(seq(0,0.5,0.025),seq(0.5,1,0.025)))
  
  # Interpolate to get rescaled values
  AUC=approx(x = aucRescale$AUC, y = aucRescale$AUCnew, xout = AUC)$y
  metric1= mean(c(COR,stability,AUC))
  metric2=mean(c(COR,stability,AUC,1-BIAS,1-MAE))
  if( metric2 > metric1) metric <- metric1 
  if (metric2 < metric1) metric <- metric2
  return(metric)
}

indexCalculation2 <- function(COR,stability,AUC,MAE,PRG,BIAS, index) {
  
  # if(AUC < 0.5) AUC<-0.5
  #  aucRescale=data.frame(AUCnew=c(rep(0,21),seq(0,1,0.05)),
  #              AUC=c(seq(0,0.5,0.025),seq(0.5,1,0.025)))
  
  # Interpolate to get rescaled values
  # AUC=approx(x = aucRescale$AUC, y = aucRescale$AUCnew, xout = AUC)$y
  metric= mean(c(COR,stability,AUC,1-MAE,1-BIAS))
  
  return(metric)
}


indexCalculation <- function(COR,stability,AUC,MAE,PRG) {
  
  if (stability < 0){stability <- 0}
  values <- c(AUC, 0, 0, 0, stability, 0, 0, 0, COR)
  metric=Lithics3D:::getTArea(values)
  metric= metric / 0.8660254
  # if(metric > 1) metric <- 1
  return(metric)
}



l$index2<-NA
lapply(1:nrow(l), function(x){
  index<-indexCalculation2(COR=l$COR[x],stability = l$stability[x],AUC=l$PRG[x],MAE=l$MAE[x], PRG=l$PRG[x], BIAS=l$BIAS[x], index=l$index[x])
  l$index2[x]<<-index
})

#l$index3<-NA
#lapply(1:nrow(l), function(x){
#  index<-indexCalculation3(COR=l$COR[x],stability = l$stability[x],AUC=l$PRG[x],MAE=l$MAE[x], PRG=l$PRG[x], BIAS=l$BIAS[x], index=l$index[x])
#  l$index3[x]<<-index
#})

#l$index<-NA
#lapply(1:nrow(l), function(x){
#  index<-indexCalculation(COR=l$COR[x],stability = l$stability[x],AUC=l$PRG[x],MAE=l$MAE[x])
#  l$index[x]<<-index
#})

if(dplyr::n_distinct(l$method)> 1) backup<-l


#l=l%>%dplyr::filter(method=="PAA")

p1<-ggplot(l, aes(x=AUC, y=trueCOR)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+facet_wrap(vars(method),nrow=1)+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

p2<-ggplot(l, aes(x=COR, y=trueCOR)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+facet_wrap(vars(method),nrow=1)+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()



p3=ggplot(l, aes(x=index2, y=trueCOR)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+ facet_wrap(vars(method),nrow=1)+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()



gridExtra::grid.arrange(p1,p2,p3,nrow=3,ncol=1)#;rm(p1,p2,p3,p4,p6)

l %>% 
  group_by(method) %>% 
  summarise(
    Rsq_index2 = R2(index2, trueCOR),
    RMSE_index2 = RMSE(index2,trueCOR),
    
    Rsq_COR = R2(COR, trueCOR),
    RMSE_COR = RMSE(COR,trueCOR),
    Rsq_AUC = R2(AUC, trueCOR),
    RMSE_AUC = RMSE(AUC,trueCOR)
  ) %>% 
  mutate_if(is.numeric, round, digits=2) 


########## fehlerbehebung


#mod=caret::train(x=l[,c("index","AUC","COR","stability","MAE","BIAS","X","indexA" )],y=l$trueCOR,method="rf",
#                trControl = trainControl())
#saveRDS(mod, "data/model.RDS")


#pred=predict(object=readRDS("data/model.RDS"),newdata=l[,c("index","AUC","COR","stability","MAE","BIAS","X","indexA" )])
#l$pred <- pred
#saveRDS(l,"data/results.RDS")

#test=l%>%dplyr::filter(noPresencePoints == 20)
