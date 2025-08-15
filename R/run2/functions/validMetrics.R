
dataPAA$validMetric<-NULL
for(i in 1:nrow(dataPAA)){
  print(i)
  df2<-dataPAA[i,]
  validMetric=c()
  if(df2$MAE > 0.35){
    validMetric=append(validMetric,1)
  }
  
  if(df2$BIAS > 0.2){
    validMetric=append(validMetric,1)
  }
  
 # if(df2$corTrainTest < 0.5){
#    validMetric=append(validMetric,1)
#  }
  
  if(abs(diff(c(df2$AUC,df2$COR))) > 0.3){
    validMetric=append(validMetric,1)
  }
  
  dataPAA$validMetric[i]<-sum(validMetric)
  rm(validMetric)
  
  
}
plot(dataPAA$diff,dataPAA$validMetric)
data$metric2<- mean(c(dataPAA$AUC,dataPAA$COR, dataPAA$stability))

dataTest=dataPAA%>%dplyr::filter(validMetric >1)

plot(dataTest$metric,dataTest$trueCor)
plot(dataPAA$metric,dataPAA$trueCor)
