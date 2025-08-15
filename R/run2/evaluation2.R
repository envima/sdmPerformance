#'@name modeling.R
#'@date 28.02.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@author 
#'@description train a model for each virtual species
#'@description code taken from the benchmark study of Valavi et al. 2023 DOI: https://doi.org/10.1111/geb.13639
#'@misc OSF repository of the study of Valavi et al. 2023: DOI: 10.17605/OSF.IO/G6DC3
#'@source https://osf.io/puk8v (of the code for the models)

# 1 - install and load packages  ####
#-----------------------------------#

library(blockCV)
library(predicts)
library(tidyverse)
library(parallel)
library(sf)

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformanceCaseStudy/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}
df=expand.grid(size=as.character(c("KNNDM","random","block1","block2","clusters")) ,
               species=c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
               points=unique(sapply(strsplit(gsub(".gpkg", "", list.files("data/run2/virtualSpeciesTrain", full.names = FALSE)), "_"), `[`, 2)),
               replicates=1:5,
               model=as.character(c("BRT")))




droplevels(df)
#vars=terra::rast("../sdmPerformance/data/variables.tif")

# create a unique name for different runs
nameRun <- paste0("run", 2)

# 4 - evaluate models ####
#------------------------#

#df=df%>%dplyr::filter(model=="RF")%>%dplyr::filter(size=="clusters")%>%dplyr::filter(testData==2)
# save test dataset or evaluate models
source("R/functions/correlationFunctions.R")


data=mclapply(1:nrow(df), function(i){
  print(i)
  if(!file.exists(paste0("data/",nameRun,"/resultsAggregated/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_points",as.character(df$points[i]),"_replicates",as.character(df$replicates[i]),".RDS"))){
    
    realDistirbution=terra::rast(paste0("../sdmPerformance/data/virtualSpecies/",as.character(df$species[i]),".tif"))
    l=list.files(paste0("data/",nameRun,"/maps/"),pattern=paste0(as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData"),full.names = T)
    l <- grep(paste0("points",as.character(df$points[i]),"_replicates",as.character(df$replicates[i])), l, value = TRUE)
    if(length(l)<2) return(NULL)
    pred=terra::rast(l);rm(l)
    pred=terra::mean(pred)
    pred=terra::mask(pred,realDistirbution)
    
    corPred=terra::layerCor(c(pred,realDistirbution),fun="cor")$correlation[[2]][1]
    
    l=list.files(paste0("data/",nameRun,"/results/"),pattern=paste0(as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData"),full.names=T)
    l <- grep(paste0("points",as.character(df$points[i]),"_replicates",as.character(df$replicates[i])), l, value = TRUE)
    
    results=do.call(rbind,lapply(l,function(x){readRDS(x)}))
    
    results=results %>%#dplyr::filter(noPresencePoints > 4)%>%
      rowwise() %>%
      mutate(
        metricSens2=mean(c(Spec,Sens,stability)),
        diff=abs(trueCor-metric),
        metricSens2_diff=abs(trueCor-metricSens2)
      ) %>%
      ungroup()
    
    nameSpecies=list.files(paste0("data/",nameRun,"/results/"),pattern=paste0(as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_testData"),full.names=F)
    nameSpecies <- grep(paste0("points",as.character(df$points[i]),"_replicates",as.character(df$replicates[i])), nameSpecies, value = TRUE)
    nameSpecies=strsplit(nameSpecies, split="_")[[1]][1]
    
    ciResults= results %>%
      group_by(method) %>%
      summarise(
        mean_metricSens2 = mean(metricSens2, na.rm = TRUE),
        mean_AUC = mean(AUC, na.rm = TRUE),
        mean_COR = mean(COR, na.rm = TRUE),
        mean_Spec = mean(Spec, na.rm = TRUE),
        mean_Sens = mean(Sens, na.rm = TRUE),
        mean_Kappa = mean(Kappa, na.rm = TRUE),
        mean_PCC = mean(PCC, na.rm = TRUE),
        mean_TSS = mean(TSS, na.rm = TRUE),
        mean_stability = mean(stability, na.rm = TRUE),
        mean_PRG = mean(PRG, na.rm = TRUE),
        mean_MAE = mean(MAE, na.rm = TRUE),
        mean_BIAS = mean(BIAS, na.rm = TRUE),
        mean_corTrainTest = mean(corTrainTest, na.rm = TRUE),
        mean_moranTest = mean(moranTest, na.rm = TRUE),
        mean_moranTrain = mean(moranTrain, na.rm = TRUE),
        mean_noPresencePoints = mean(noPresencePoints, na.rm = TRUE),
        mean_metric = mean(metric, na.rm = TRUE),
        
        lower_ci = if (sum(!is.na(metricSens2)) >= 2 && sd(metricSens2, na.rm = TRUE) > 0) {
          t.test(metricSens2)$conf.int[1]
        } else {
          NA_real_
        },
        upper_ci = if (sum(!is.na(metricSens2)) >= 2 && sd(metricSens2, na.rm = TRUE) > 0) {
          t.test(metricSens2)$conf.int[2]
        } else {
          NA_real_
        },
        species=nameSpecies,
        points=first(points),
        replicates=first(replicate),
        model=first(model),
        size = first(size),
        testData = first(testData),
        trueCor=first(corPred)
      )
    if(!dir.exists(paste0("data/",nameRun,"/resultsAggregated"))) dir.create(paste0("data/",nameRun,"/resultsAggregated"),recursive=T)
    saveRDS(ciResults, paste0("data/",nameRun,"/resultsAggregated/",as.character(df$species[i]),"_",as.character(df$size[i]),"_",as.character(df$model[i]),"_points",as.character(df$points[i]),"_replicates",as.character(df$replicates[i]),".RDS"))
    return(ciResults)
  }
},mc.cores=nCores)

data=do.call(rbind,data)
saveRDS(data, "data/run1/results_ci_metricSens2.RDS")
data=readRDS("data/run1/results_ci_metricSens2.RDS")




data %>%
  mutate(order = rank(trueCor)) %>%  # keep it numeric
  ggplot(aes(x = order, y = trueCor)) +
  geom_point() +  #ylim(0.5,1.5)+
  geom_point(aes(y = mean_metricSens2, colour = "metric"), shape = 2,color="blue") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "lightblue", alpha = 0.4) +
  facet_wrap(vars(method), scales = "free_x") +
  labs(x = "Ordered by trueCor", y = "trueCor") 




# ###################################################
#
#
# Show the probem:
#
#
######################################################


results=list.files(paste0("data/",nameRun,"/results1/"),full.names=T)
results=do.call(rbind,lapply(results,function(x){readRDS(x)}))

results%>%dplyr::filter(method %in% c("PA","PBG"))%>%
  ggplot(aes(y=AUC,x=size))+
  geom_boxplot()+facet_wrap(vars(method,model))+xlab("Method used to separate the folds")


results%>%dplyr::filter(method %in% c("PA","PBG"))%>%
  ggplot(aes(y=COR,x=size))+
  geom_boxplot()+facet_wrap(vars(method,model))+xlab("Method used to separate the folds")

#####################


library(ggplot2)
library(dplyr)
library(ggpmisc)

results$identifier <- paste(results$model, results$size,results$testData,sep="_")

# do for
for (i in c("AUC","MAE","Spec","Sens")){}

results %>%
  dplyr::filter(method %in% c("PA", "PBG")) %>%
  ggplot(aes(x = AUC, y = trueCor)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue",fill="lightblue") +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    size = 4
  ) +
  facet_wrap(vars(method))


results %>%
  dplyr::filter(method %in% c("PA")) %>%
  ggplot(aes(x = AUC, y = trueCor)) +
  geom_point() +xlim(0.5,1.2)+
  geom_smooth(method = "lm", se = T, color = "blue",fill="lightblue") +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    size = 4
  ) +
  #facet_wrap(vars(method))+
  geom_text(aes(label=ifelse(COR>0.8,as.character(identifier),'')),hjust=0,vjust=0)

# get the "best model"

bestModel=results%>%dplyr::filter(method=="PA")%>%dplyr::filter(AUC==1)
bestModel=bestModel%>%dplyr::filter(MAE==min(bestModel[bestModel$method=="PA",]$MAE))

standard=results%>%dplyr::filter(method=="PA")%>%dplyr::filter(model=="Maxent")%>%dplyr::filter(size=="block1")%>%
  dplyr::filter(testData==1)


# ##################################################################


# What do we do about this?



####################################################################

results2=results %>%#dplyr::filter(noPresencePoints > 4)%>%
  rowwise()%>%
  mutate(
    metricSens2=mean(c(Sens,Spec,stability)))

results2 %>%
  dplyr::filter(method %in% c("PA", "PBG")) %>%
  #dplyr::filter(noPresencePoints > 4)%>%
  ggplot(aes(x = metricSens2, y = trueCor)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue",fill="lightblue") +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    size = 4
  ) +
  facet_wrap(vars(method,noPresencePoints))

bestModel=results2%>%dplyr::filter(metricSens2==max(results2$metricSens2))
worstModel=results2%>%dplyr::filter(metricSens2==min(results2$metricSens2))



