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

setwd("M:/user/bald/SDM/sdmPerformance/")

df=expand.grid(vs=gsub(".tif","",list.files("data/paRaster/", pattern=".tif", full.names=F)) ,
               points=c(40,120,200,400),
               replicates=1,
               cv="random",
               testData=1:6,
               method=NA,
               index=NA,
               AUC=NA,
               COR=NA,
               stability=NA)


#df=df%>%dplyr::filter(df$vs == "VS01")

if(!dir.exists("data/index2/")) dir.create("data/index2")

#mclapply(1:nrow(df), function(s){
lapply(1:nrow(df), function(s){
  print(s)
  if(!file.exists(sprintf("data/index2/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))){
    #for (s in 1:nrow(df)){
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
                                    noPointsTesting=NA)
    dfPAA=df2
    dfPAA$index<- index$indexPAA$metric
    dfPAA$AUC<- index$indexPAA$AUC
    dfPAA$COR<- index$indexPAA$COR
    dfPAA$PRG<- index$indexPAA$PRG
    dfPAA$MAE<- index$indexPAA$MAE
    dfPAA$BIAS<- index$indexPAA$BIAS
    dfPAA$stability<- index$indexPAA$stability
    dfPAA$method <- "PAA"
    
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
    }
    dfPBG=df2
    dfPBG$index<- index$indexPBG$metric
    dfPBG$AUC<- index$indexPBG$AUC
    dfPBG$COR<- index$indexPBG$COR
    dfPBG$stability<- index$indexPBG$stability
    dfPBG$PRG<- index$indexPBG$PRG
    dfPBG$MAE<- index$indexPBG$MAE
    dfPBG$BIAS<- index$indexPBG$BIAS
    dfPBG$method <- "PBG"
    
    if (is.na(index$indexPA)[1]){
      df2=rbind( dfPAA, dfPBG)
    } else {df2=rbind(dfPA, dfPAA, dfPBG)}
    
    rm(index, p,vs, dfPA, dfPAA, dfPBG)
    saveRDS(df2, sprintf("data/index2/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))
  }
})
#}, mc.cores=45)


l=list.files("data/results_1000/", full.names=T)
l=do.call(rbind, lapply(l, function(x){
  data=readRDS(x)
  return(data)}))

#data=l%>%dplyr::mutate(diffPA=abs(COR-indexPA),
#                       diffPAA=abs(COR-indexPAA),
#                       diffPBG=abs(COR-indexPBG))



#dataPA=data.frame(value=data$diffPA, method="indexPA")
#dataPAA=data.frame(value=data$diffPAA, method="indexPAA")
#dataPBG=data.frame(value=data$diffPBG, method="indexPBG")
#dataPlots=rbind(dataPA, dataPAA, dataPBG);rm(dataPA,dataPAA,dataPBG)

#data=l %>% dplyr::filter(method=="PAA")



#long <- data %>% 
#  pivot_longer(
#    cols = c(  "index" ,     "AUC" ,       "COR"    ,"stability",  "trueCOR" ), 
#    names_to = "metric",
#    values_to = "value"
#  )

# Plot
#long %>% #dplyr::filter(method==PA)%>%
#  ggplot( aes(x=metric,y=value, color=metric, shape=metric)) + 
#  geom_point(size=4) +facet_wrap(vars(vs,points))


#x=l#%>%dplyr::filter(method=="PBG")#%>%dplyr::filter(points==40)
#cor(x$COR, x$trueCOR)
#cor(x$AUC, x$trueCOR)
#cor(x$index, x$trueCOR)
#cor(x$stability , x$COR)

#x=l%>%dplyr::filter(method=="PA")
#plot(x$index, x$trueCOR)



l[l < 0] <- 0
#l=l%>%dplyr::filter(method=="PAA")

# Mean absoulte error of outliers
#observed=l$trueCOR
#predicted=l$COR
#residuals <- observed - predicted
#threshold <- 2 * sd(residuals)
#outliers <- abs(residuals) > threshold
#mae_outliers <- mean(abs(residuals[outliers]))


indexCalculation2 <- function(COR,stability,AUC,MAE,PRG) {
  

 metric=mean(c(COR,AUC,stability))
 
  
  return(metric)
}

indexCalculation <- function(COR,stability,AUC,MAE) {
 
  #if (stability < 0){stability <- 0}
  values <- c(COR, 0, 0, 0, stability, 0, 0, 0, AUC)
  metric=Lithics3D:::getTArea(values)
  metric= metric / 0.8660254
 # if(metric > 1) metric <- 1
  return(metric)
}

l$index3<-NA
lapply(1:nrow(l), function(x){
  index<-indexCalculation2(COR=l$COR[x],stability = l$stability[x],AUC=l$PRG[x],MAE=l$MAE[x], PRG=l$PRG[x])
  l$index3[x]<<-index
})

l$index<-NA
lapply(1:nrow(l), function(x){
  index<-indexCalculation(COR=l$COR[x],stability = l$stability[x],AUC=l$PRG[x],MAE=l$MAE[x])
  l$index[x]<<-index
})

p1<-ggplot(l, aes(x=index, y=trueCOR,color=points)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

p2<-ggplot(l, aes(x=index3, y=trueCOR,color=points)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

p3=ggplot(l, aes(x=COR, y=trueCOR,color=points))+
  geom_point()+ylim(0,1)+xlim(0,1)+facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()


p4=ggplot(l, aes(x=AUC, y=trueCOR,color=points)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+ facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

p5<-ggplot(l, aes(x=1-MAE, y=trueCOR,color=points)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+ facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

gridExtra::grid.arrange(p2,p3,p4,nrow=3,ncol=1)

#x=l%>%dplyr::filter(method=="PAA")

summ <- l %>% 
  group_by(method) %>% 
  summarise(Rsq_index = R2(index, trueCOR),
            RMSE_index = RMSE(index,trueCOR),
            Rsq_index3 = R2(index3, trueCOR),
            RMSE_index3 = RMSE(index3,trueCOR),
            Rsq_COR = R2(COR, trueCOR),
            RMSE_COR = RMSE(COR,trueCOR),
            Rsq_AUC = R2(AUC, trueCOR),
            RMSE_AUC = RMSE(AUC,trueCOR)
  ) %>% 
  mutate_if(is.numeric, round, digits=2) 
summ

#plot(log(seq(0,1,0.1)),seq(0,1,0.1) )

