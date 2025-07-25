
meanBestModelFunction <- function(s, # species name
                                  p, # number of points
                                  rep, # number of replicates
                                  m, # method used for calculation of metrics (PA, PAA, PBG)
                                  evaluationIndex # e.g. AUC, metric, COR, TSS...
){
  if(!file.exists(paste0("data/run2/resultsAggregated/",m,"_",evaluationIndex,"/", s, "_",p,"_", rep, ".RDS"))){
    #load result file
    results=readRDS("data/run2/results.RDS")%>%
      dplyr::filter(method==m)%>%
      dplyr::filter(species==s, points==p, replicate==rep)
    
    
    results=results %>%
      dplyr::mutate(
        metric = rowMeans(cbind(Spec, COR, PCC, 1 - MAE), na.rm = TRUE)
      )
    
    
    # get best 90% of models
    bestModels = results%>%dplyr::filter(.data[[evaluationIndex]] >= quantile(results[[evaluationIndex]], probs=c(0.9)))
    
    
    predictions=lapply(1:nrow(bestModels), function(x){
      pred=terra::rast(sprintf("data/run2/maps/%s_%s_%s_testData%s_points%s_replicates%s.tif", bestModels$species[x],
                               bestModels$size[x],
                               bestModels$model[x],
                               bestModels$testData[x],
                               bestModels$points[x],
                               bestModels$replicate[x]))
      return(pred)
    })
    
    
    predMean=terra::mean(terra::rast(predictions), na.rm=T)
    realDistirbution=terra::rast(paste0("data/virtualSpecies/",s,".tif"))
    predMean=terra::mask(predMean, realDistirbution)
    x=terra::rast(list(predMean, realDistirbution))
    
    
    
    
    
    bestModelMean <- bestModels%>%dplyr::summarise(metricMean = mean(metric,na.rm=T),
                                                   AUCMean = mean(AUC,na.rm=T),
                                                   CORMean = mean(COR,na.rm=T),
                                                   SpecMean = mean(Spec,na.rm=T),
                                                   SensMean = mean(Sens,na.rm=T),
                                                   KappaMean = mean(Kappa,na.rm=T),
                                                   PCCMean = mean(PCC,na.rm=T),
                                                   TSSMean = mean(TSS,na.rm=T),
                                                   PRGMean = mean(PRG,na.rm=T),
                                                   MAEMean = mean(MAE,na.rm=T),
                                                   BIASMean = mean(BIAS,na.rm=T),
                                                   trueCor_Mean=mean(trueCor, na.rm=T)
    )%>%dplyr::mutate(stability=stabilityRasters(r=terra::rast(predictions)))
    
    
    
    index=mean(c(bestModelMean$SpecMean, bestModelMean$CORMean, bestModelMean$PCCMean,1-  bestModelMean$MAEMean,bestModelMean$stability),na.rm=T)
    bestModelMean$metricMean<- index
    
    
    df2=data.frame(species=s,
                   points=p,
                   replicates=rep,
                   trueCorMean = terra::layerCor(x, fun="cor")$correlation[2])
    
    df2[[evaluationIndex]]= bestModelMean[[paste0(evaluationIndex,"Mean")]]
    if(!dir.exists(paste0("data/run2/resultsAggregated/",m,"_",evaluationIndex))) dir.create(paste0("data/run2/resultsAggregated/",m,"_",evaluationIndex), recursive=T)
    saveRDS(df2,paste0("data/run2/resultsAggregated/",m,"_",evaluationIndex,"/", s, "_",p,"_", rep, ".RDS"))
    
    return(df2)
  } # end file exists
}