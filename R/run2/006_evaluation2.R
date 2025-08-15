#'@name 006_evaluation2.R
#'@date 24.07.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]

#'@description 
#'@description 


# 1 - install and load packages  ####
#-----------------------------------#

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=40
}

source("R/functions/performanceStabilityIndexReplicates.R")
source("R/functions/meanBestModelFunction.R")
library(tidyverse)
library(parallel)
library(ggplot2)


set.seed(12345)

# 2 - all model combinations ####
#-------------------------------#

df=expand.grid(species=c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06","VS07", "VS08", "VS09", "VS10"),
               points=unique(sapply(strsplit(gsub(".gpkg", "", list.files("data/run2/virtualSpeciesTrain", full.names = FALSE, pattern=".gpkg")), "_"), `[`, 2)),
               replicates=1:5)


# 3 - start calculation ####
#--------------------------#


data=mclapply(1:nrow(df), function(i){
  print(i)
  df2=meanBestModelFunction(s=as.character(df$species[i]),
                            p=as.character(df$points[i]),
                            rep=as.character(df$replicate[i]),
                            m="PBG", 
                            evaluationIndex="metric")
},mc.cores=nCores)

data=mclapply(1:nrow(df), function(i){
  print(i)
  df2=meanBestModelFunction(s=as.character(df$species[i]),
                            p=as.character(df$points[i]),
                            rep=as.character(df$replicate[i]),
                            m="PBG", 
                            evaluationIndex="AUC")
},mc.cores=nCores)

data=mclapply(1:nrow(df), function(i){
  print(i)
  df2=meanBestModelFunction(s=as.character(df$species[i]),
                            p=as.character(df$points[i]),
                            rep=as.character(df$replicate[i]),
                            m="PAA", 
                            evaluationIndex="metric")
},mc.cores=nCores)

data=mclapply(1:nrow(df), function(i){
  print(i)
  df2=meanBestModelFunction(s=as.character(df$species[i]),
                            p=as.character(df$points[i]),
                            rep=as.character(df$replicate[i]),
                            m="PA", 
                            evaluationIndex="metric")
},mc.cores=nCores)

# 4 - read and plot data ####
#---------------------------#


data=list.files("data/run2/resultsAggregated/PAA_metric", pattern=".RDS", full.names=T)
data=do.call( "rbind",lapply(data, readRDS))


data%>%ggplot(aes(x=metric, y=trueCorMean))+
  geom_point()+xlim(0,1)+ylim(0,1)



