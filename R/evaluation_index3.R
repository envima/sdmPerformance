#'@name evaluation.R
#'@date 27.05.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@description train a model for each virtual species
#'@description code taken from the benchmark study of Valavi et al. 2023 DOI: https://doi.org/10.1111/geb.13639
#'@misc OSF repository of the study of Valavi et al. 2023: DOI: 10.17605/OSF.IO/G6DC3
#'@source https://osf.io/puk8v (of the code for the models)



source("R/functions/performanceStabilityIndexReplicates.R")
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
               testData=1:6)


# create output dir
if(!dir.exists("data/output/Maxent/resultsReplicates")) dir.create("data/output/Maxent/resultsReplicates")

mclapply(1:nrow(df), function(s){
  print(s)
  if(!file.exists(sprintf("data/output/Maxent/resultsReplicates/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))){
    df2=df[s,]
    # calculate true correlation
    p=list.files("data/output/Maxent/predictions/", pattern=sprintf("Maxent_%s_%s_%s_%s_testData%s_cvFold",df$vs[s], df$points[s],df$replicates[s],df$cv[s],df$testData[s]),full.names=T)
    p=terra::rast(p)
    vs=terra::rast(paste0("data/virtualSpecies/",df$vs[s],".tif"))
    df2$trueCOR<- terra::layerCor(c(vs, terra::mean(p)),fun="cor")$correlation[1,2]
    if (df2$trueCOR < 0) df2$trueCOR <-0
    # calculate performance metrics
    index=performanceStabilityIndex(x=p,
                                    absence = sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==0)%>%dplyr::filter(fold == df$testData[s]),
                                    presence=sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==1)%>%dplyr::filter(fold == df$testData[s]),
                                    background=TRUE,
                                    trainingData=sf::read_sf(sprintf("data/processed/dataPartition_random/PA/%s_%s_%s.gpkg",df$vs[s],df$points[s],df$replicates[s]))%>%dplyr::filter(Real==1)%>%dplyr::filter(fold != df$testData[s]),
                                    aa=TRUE,
                                    environmentalVariables=terra::rast("data/variables.tif"),
                                    noPointsTesting=NA)
    # bring everything together in one df
    index$vs <- df2$vs
    index$points <- df2$points
    index$replicates <- df2$replicates
    index$cv <- df2$cv
    index$testData <- df2$testData
    index$trueCOR <- df2$trueCOR
    index$method <- rownames(index)
    
    # save results
    saveRDS(index, sprintf("data/output/Maxent/resultsReplicates/Maxent_%s_%s_%s_%s_testData%s.RDS",df$vs[s], df$points[s],df$replicates[s],df$cv[2],df$testData[s]))
    rm(index, p,vs,df2);gc()
  }
  
}, mc.cores=45)


if (!file.exists("data/output/Maxent/resultsReplicates.RDS")){
  l=list.files("data/output/Maxent/resultsReplicates", full.names=T,pattern=".RDS")
  l=do.call(rbind, lapply(l, function(x){
    data=readRDS(x)
    return(data)}))
  saveRDS(l, "data/output/Maxent/results.RDS")
};rm(l)



# 3 - Plot results ####
#---------------------#








data=readRDS("data/output/Maxent/results.RDS")
data$method=gsub("index","",data$method)


model_AUC <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ AUC, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["AUC"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_COR <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ COR, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["COR"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_metric <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ metric, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["metric"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

p1<- ggplot(data, aes(x = AUC, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_AUC,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")


p2<- ggplot(data, aes(x = COR, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_COR,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")

p3<- ggplot(data, aes(x = metric, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_metric,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")



gridExtra::grid.arrange(p1,p2,p3,nrow=3,ncol=1)#;rm(p1,p2,p3,p4,p6)

data %>% 
  group_by(method) %>% 
  summarise(
    Rsq_index2 = R2(metric, trueCOR),
    RMSE_index2 = RMSE(metric,trueCOR),
    
    Rsq_COR = R2(COR, trueCOR),
    RMSE_COR = RMSE(COR,trueCOR),
    Rsq_AUC = R2(AUC, trueCOR),
    RMSE_AUC = RMSE(AUC,trueCOR)
  ) %>% 
  mutate_if(is.numeric, round, digits=2) 


