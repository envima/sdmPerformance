#'@name plots.R
#'@date 21.07.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@description plots the results

# 1 - install and load packages  ####
#-----------------------------------#

library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(parallel)
library(hexbin)



if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformanceCaseStudy/")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

# 2 - load data ####
#------------------#

data=readRDS("data/run2/results.RDS")

data <- data %>% filter(is.finite(metric))%>%
  mutate(trueCor = ifelse(trueCor < 0, 0, trueCor))


# 3 - main plot ####
#------------------#

# Calculate R², RMSE, and slope per method
metrics_df <- data %>%
 # filter(noPresencePoints > 5) %>%
  group_by(method) %>%
  summarise(
    r_squared = summary(lm(trueCor ~ metric))$r.squared,
    rmse = sqrt(mean((trueCor - predict(lm(trueCor ~ metric)))^2, na.rm = TRUE)),
    slope = coef(lm(trueCor ~ metric))[["metric"]],
    xpos = 0,
    ypos = 1
  )

# Plot with text annotations for R² and RMSE
data %>%
  #filter(noPresencePoints > 5) %>%
  ggplot(aes(x = metric, y = trueCor)) +
  # geom_point(size = 0.005) +
  geom_hex(bins = 140, alpha=1) +
  scale_fill_continuous(type = "viridis", name = "Point density") +
  geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
  facet_wrap(vars(method)) +
  geom_text(
    data = metrics_df,
    aes(x = xpos, y = ypos, 
        label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope,2))),
    inherit.aes = FALSE,
    size = 5,
    hjust = 0
  )+xlim(0,1)+ylim(-0.01,1)+
  xlab("Metric")+ylab("Pearson correlation between suitability \nraster and predcition map of virtual species")+
  theme_minimal(base_size = 15) #+
#  theme(
#    legend.title = element_text(size = 14, face = "bold")
#  )

ggplot2::ggsave(filename="")

#######################################################################

# Define the metric names
metric_names <- c("metric", "AUC", "COR", "Spec", "Sens", "Kappa", 
                  "PCC", "TSS", "PRG", "MAE", "BIAS")

# Create empty list to store plots
plot_list <- list()

# Loop over each metric name
for (m in metric_names) {
  # Filter and prepare data
  df <- data %>%
    filter(noPresencePoints > 5, is.finite(.data[[m]]), is.finite(trueCor))
  
  # Calculate RMSE and R² for each method
  metrics_df <- df %>%
    group_by(method) %>%
    summarise(
      r_squared = summary(lm(trueCor ~ .data[[m]]))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(trueCor ~ .data[[m]])))^2, na.rm = TRUE)),
      xpos = quantile(.data[[m]], 0.05, na.rm = TRUE),
      ypos = quantile(trueCor, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create ggplot
  p <- ggplot(df, aes(x = .data[[m]], y = trueCor)) +
    geom_point(size = 0.1) +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(method)) +
    geom_text(data = metrics_df,
              aes(x = xpos, y = ypos,
                  label = paste0("R² = ", round(r_squared, 2), 
                                 "\nRMSE = ", round(rmse, 2))),
              inherit.aes = FALSE,
              size = 3, hjust = 0) +
    labs(title = paste("trueCor vs", m), x = m, y = "trueCor")
  
  # Save to list
  plot_list[[m]] <- p
}














#------------------------------------------------------










data$ID <- paste(data$species, data$model, data$method, data$size, data$points, data$replicates, sep="_")
dataCluster = data%>%dplyr::filter(size=="clusters")



data %>%
  dplyr::slice_sample(n=100)%>%
  dplyr::filter(method == "PAA") %>%
  # dplyr::arrange(trueCor) %>%
  #  dplyr::mutate(ID_ordered = factor(ID, levels = ID)) %>%
  # ggplot(aes(x = ID_ordered)) +
  ggplot(aes(x = ID)) +
  # facet_wrap(vars(size))+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "grey80", alpha = 0.4) +
  geom_point(aes(y = trueCor), color = "steelblue") +
  geom_point(aes(y = mean_metricSens2), color = "darkred") +
  labs(x = "Modell-ID", y = "Wert", title = "trueCor (blau) und mean_metric (rot) mit CI-Band") +
  theme_minimal()


geom_line(aes( y = bestAUC, group = 1, color = "bestAUC"))+
  geom_line(aes( y = bestCOR, group = 1, color = "bestCOR"))+
  geom_line(aes( y = bestMetric, group = 1, color = "bestMetric"))+
  #geom_point() +
  labs(x = "Model Identifier", y = "model rank") #+
#  theme_minimal() +
#coord_flip()








data$diff <- abs(data$COR-data$trueCor)

#data=data%>%dplyr::filter(method=="indexPA")

#
#   How good are the models with different spatial cv ? 
#

# filtering should not matter as each true cor value is in there three times. just choose one.
#data%>%dplyr::filter(method=="PAA")%>%
#ggplot( aes(x = reorder(size, trueCor, FUN = median), y = trueCor)) + 
#  geom_boxplot() #+  geom_jitter()




dataPAA <- data %>%dplyr::filter(noPresencePoints >= 5)%>%
  filter(method == "PAA") %>%
  mutate(
    bestModel = rank(-trueCor, ties.method = "min"),
    bestAUC = rank(-AUC, ties.method = "min"),
    bestCOR = rank(-COR, ties.method = "min"),
    bestMetric = rank(-metric, ties.method = "min"),
    modelID = bestModel
  ) %>%
  rowwise() %>%
  mutate(
    maeMean = mean(c(BIAS, MAE, 1 - stability), na.rm = TRUE),
    maeIntervallMin = metric - maeMean,
    maeIntervallMax = metric + maeMean,
    biasRange = abs(maeIntervallMin - maeIntervallMax),
    newMetric=metric+mean(c(BIAS, MAE, 1 - stability)/3, na.rm = TRUE),
    metric2=median(c(AUC,COR,stability,PRG)),
    metric3=mean(c(AUC,COR,stability,TSS)),
    metricSens=mean(c(Spec,Sens)),
    metricSens2=mean(c(Spec,Sens,stability)),
    newMetric6=metric+mean(c(BIAS, MAE, 1 - stability)/2, na.rm = TRUE),
    newMetric2=metric-mean(c(BIAS, MAE, 1 - stability)/3, na.rm = TRUE),
    newMetric3=metric+mean(c(BIAS, MAE, 1 - stability)/2, na.rm = TRUE),
    newMetric4=metric+mean(c(BIAS, MAE, 1 - stability)/4, na.rm = TRUE),
    newMetric5=mean(c(metric+mean(c(BIAS, MAE, 1 - stability )/3, na.rm = TRUE),PRG)),
    newMetric7=metricSens2-abs(log(BIAS)/100),
    hibridMetric=(0.4*AUC + 0.3*COR + 0.3*Kappa),
    
    diff=abs(trueCor-metricSens2)
  ) %>%
  ungroup()
rownames(dataPAA)<-1:nrow(dataPAA)


#plot(dataPAA$trueCor, dataPAA$corTrainTest)
#plot(dataPAA$trueCor, dataPAA$moranTrain)
#plot(dataPAA$trueCor, dataPAA$moranTest)



dataPAA %>%
  #dplyr::filter(validMetric < 1)%>%
  dplyr::filter(noPresencePoints > 4)%>%
  # dplyr::filter(corTrainTest > 0.4)%>%
  #  dplyr::filter(MAE < 0.35)%>%
  #dplyr::filter(BIAS<0.15)%>%
  ggplot(aes(x = reorder(modelID, trueCor), y = trueCor, group = 1)) +#facet_wrap(vars(noPresencePoints))+
  geom_line() +#ylim(.55,1)+
  # geom_line(aes( y = AUC, group = 1, colour = "AUC"))+
  #  geom_line(aes( y = MAE, group = 1, colour = "MAE"))+
  # geom_line(aes( y = BIAS, group = 1, colour = "BIAS"))+
  #  geom_line(aes( y = metricSens, group = 1, colour = "metricSens"))+
  geom_line(aes( y = metricSens2, group = 1, colour = "metricSens2"))+
  geom_line(aes( y = moranTest, group = 1, colour = "corTrainTest"))+
  #  geom_line(aes( y = metric, group = 1, colour = "metric"))+
  # geom_line(aes( y = COR, group = 1, colour = "COR"))+
  #geom_line(aes( y = newMetric7, group = 1, colour = "newMetric7"))+
  # geom_point(aes( y = newMetric3, group = 1, colour = "newMetric3"))+
  # geom_point(aes( y = newMetric6, group = 1, colour = "newMetric6"))+
  #geom_line(aes( y = newMetric5, group = 1, colour = "newMetric5"))+
  # geom_ribbon(aes(ymin = maeIntervallMin, ymax = maeIntervallMax), fill = "lightblue", alpha = 0.4) +
  # geom_line(aes(y = trueCor), color = "darkblue", size = 1) +
  labs(x = "Model Identifier", y = "Value", title = "True Correlation and MAE Interval") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


plot(dataPAA$diff,log(dataPAA$BIAS))
cor(dataPAA$diff,log(dataPAA$BIAS))
plot(dataPAA$diff,dataPAA$MAE)
cor(dataPAA$diff,dataPAA$MAE)
plot(dataPAA$diff,log(dataPAA$moranTest))
cor(dataPAA$diff,log(dataPAA$moranTest))
plot(dataPAA$diff,dataPAA$moranTrain)
cor(dataPAA$diff,dataPAA$moranTrain)
plot(dataPAA$diff,dataPAA$corTrainTest)
cor(dataPAA$diff,dataPAA$corTrainTest)
plot(dataPAA$COR,dataPAA$diff )
cor(dataPAA$trueCor,dataPAA$COR)



dataTest=dataPAA%>%dplyr::filter(noPresencePoints >= 5)%>%dplyr::filter(size=="KNNDM")%>%dplyr::filter(testData==2)#dplyr::filter(MAE > 0.4)
dataTest$metric_corrected<-dataTest$metric+0.095
dataTest[c("trueCor","metric_corrected","metric")]
##############################################################################





####################################

dataPAA %>%dplyr::filter(noPresencePoints > 5)%>%
  ggplot(aes(x = reorder(modelID, bestModel), y = bestModel, group = 1)) +#facet_wrap(vars(noPresencePoints))+
  geom_line() +
  geom_line(aes( y = bestAUC, group = 1, color = "bestAUC"))+
  geom_line(aes( y = bestCOR, group = 1, color = "bestCOR"))+
  geom_line(aes( y = bestMetric, group = 1, color = "bestMetric"))+
  #geom_point() +
  labs(x = "Model Identifier", y = "model rank") #+
#  theme_minimal() +
#coord_flip()

dataPAA %>%dplyr::filter(noPresencePoints > 5)%>%
  ggplot(aes(x = reorder(modelID, trueCor), y = trueCor, group = 1)) +#facet_wrap(vars(noPresencePoints))+
  geom_line() +
  geom_line(aes( y = AUC, group = 1, colour = "AUC"))+
  geom_line(aes( y = COR, group = 1, colour = "COR"))+
  geom_line(aes( y = metric, group = 1, colour = "metric"))+
  #geom_point() +
  labs(x = "Model Identifier", y = "metric and true COR")# +
#  theme_minimal()# +
# coord_flip()

dataPAA %>%#dplyr::filter(noPresencePoints > 20)%>%
  ggplot(aes(x = metric, y = trueCor, group = 1)) +facet_wrap(vars(noPresencePoints))+
  geom_point()

dataMAE <- dataPAA%>%dplyr::filter(MAE < 0.3)

library(viridis)
dataMAE %>%#dplyr::filter(noPresencePoints > 10)%>%
  ggplot(aes(x = metric, y = trueCor, colour = noPresencePoints)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  scale_color_viridis_c()  # Use _d() if the variable is categorical




test=dataPAA%>%dplyr::group_by(model,size)%>%dplyr::summarise()

test %>%#dplyr::filter(noPresencePoints > 20)%>%
  ggplot(aes(x = reorder(modelID, trueCor), y = trueCor, group = 1)) +#facet_wrap(vars(noPresencePoints))+
  geom_line() +
  #geom_line(aes(y=noPresencePoints))+
  geom_line(aes( y = AUC, group = 1, colour = "AUC"))+
  geom_line(aes( y = COR, group = 1, colour = "COR"))+
  geom_line(aes( y = metric, group = 1, colour = "metric"))+
  #geom_point() +
  labs(x = "Model Identifier", y = "metric and true COR")# +
#  theme_minimal()# +
# coord_flip()

# this shows the real model performance
data%>%dplyr::filter(method=="PA")%>%
  ggplot( aes(y=reorder(identifier,trueCor), x=trueCor)) + ylab("Model specifications")+
  geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


data%>%dplyr::filter(method=="PA")%>%
  ggplot( aes(y=reorder(identifier,AUC), x=AUC)) + ylab("AUC on Presence-Absence data")+ xlab("Model specifications") +
  geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

data%>%dplyr::filter(method=="PBG")%>%
  ggplot( aes(y=reorder(identifier,AUC), x=AUC)) + xlab("Model specifications")+ ylab("AUC on Presence-Background data")+ 
  geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

data%>%dplyr::filter(method %in% c("PA"))%>%
  ggplot( aes(y=reorder(identifier,COR), x=COR)) + ylab("Model specifications")+ xlab("Pearson correlation on Presence-Absence data")+ 
  geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



