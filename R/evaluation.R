source("R/functions/performanceStabilityIndex.R")
library(tidyverse)
library(parallel)
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

if(!dir.exists("data/results_1000/")) dir.create("data/results_1000")

mclapply(1:nrow(df), function(s){
  print(s)
  if(!file.exists(sprintf("data/results_1000/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))){
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
                                    noPointsTesting=1000)
    dfPAA=df2
    dfPAA$index<- index$indexPAA$metric
    dfPAA$AUC<- index$indexPAA$AUC
    dfPAA$COR<- index$indexPAA$COR
    dfPAA$stability<- index$indexPAA$stability
    dfPAA$method <- "PAA"
    
    if (!is.na(index$indexPA)[1]){
      dfPA=df2
      dfPA$index<- index$indexPA$metric
      dfPA$AUC<- index$indexPA$AUC
      dfPA$COR<- index$indexPA$COR
      dfPA$stability<- index$indexPA$stability
      dfPA$method <- "PA"
    }
    dfPBG=df2
    dfPBG$index<- index$indexPBG$metric
    dfPBG$AUC<- index$indexPBG$AUC
    dfPBG$COR<- index$indexPBG$COR
    dfPBG$stability<- index$indexPBG$stability
    dfPBG$method <- "PBG"
    
    if (is.na(index$indexPA)[1]){
      df2=rbind( dfPAA, dfPBG)
    } else {df2=rbind(dfPA, dfPAA, dfPBG)}
    
    rm(index, p,vs, dfPA, dfPAA, dfPBG)
    saveRDS(df2, sprintf("data/results_1000/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))
  }
  #}
}, mc.cores=45)


l=list.files("data/results/", full.names=T)
l=do.call(rbind, lapply(l, function(x){
  data=readRDS(x)
  return(data)}))

data=l%>%dplyr::mutate(diffPA=abs(COR-indexPA),
                       diffPAA=abs(COR-indexPAA),
                       diffPBG=abs(COR-indexPBG))

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

dataPA=data.frame(value=data$diffPA, method="indexPA")
dataPAA=data.frame(value=data$diffPAA, method="indexPAA")
dataPBG=data.frame(value=data$diffPBG, method="indexPBG")
dataPlots=rbind(dataPA, dataPAA, dataPBG);rm(dataPA,dataPAA,dataPBG)

data=l %>% dplyr::filter(method=="PAA")

library(tidyr)

long <- data %>% 
  pivot_longer(
    cols = c(  "index" ,     "AUC" ,       "COR"    ,"stability",  "trueCOR" ), 
    names_to = "metric",
    values_to = "value"
  )

# Plot
long %>% #dplyr::filter(method==PA)%>%
  ggplot( aes(x=metric,y=value, color=metric, shape=metric)) + 
  geom_point(size=4) +facet_wrap(vars(vs,points))


x=l#%>%dplyr::filter(method=="PBG")#%>%dplyr::filter(points==40)
cor(x$COR, x$trueCOR)
cor(x$AUC, x$trueCOR)
cor(x$index, x$trueCOR)
cor(x$stability , x$COR)

x=l%>%dplyr::filter(method=="PA")
plot(x$index, x$trueCOR)

# Library
library(ggplot2)
library(hrbrthemes)
library(caret)
library(tidyverse)
library(ggpmisc)

l[l < 0] <- 0


summ <- l %>% 
  group_by(method) %>% 
  summarise(Rsq_index = R2(index, trueCOR),
            RMSE_index = RMSE(index,trueCOR),
            Rsq_COR = R2(COR, trueCOR),
            RMSE_COR = RMSE(COR,trueCOR),
            Rsq_AUC = R2(AUC, trueCOR),
            RMSE_AUC = RMSE(AUC,trueCOR),
            MAE_COR =caret::MAE(COR, trueCOR),
            MAE_index =caret::MAE(index, trueCOR),
            MAE_AUC =caret::MAE(AUC, trueCOR)
            ) %>% 
  mutate_if(is.numeric, round, digits=2) 

# Mean absoulte error of outliers
observed=l$trueCOR
predicted=l$COR
residuals <- observed - predicted
threshold <- 2 * sd(residuals)
outliers <- abs(residuals) > threshold
mae_outliers <- mean(abs(residuals[outliers]))



p1<-ggplot(l, aes(x=index, y=trueCOR)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

p2=ggplot(l, aes(x=COR, y=trueCOR))+
  geom_point()+ylim(0,1)+xlim(0,1)+facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()


p3=ggplot(l, aes(x=AUC, y=trueCOR)) + 
  geom_point() +ylim(0,1)+xlim(0,1)+ facet_wrap(vars(method))+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()


p=gridExtra::grid.arrange(p1,p2,p3,nrow=3,ncol=1)
