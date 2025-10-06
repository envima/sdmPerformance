#'@name plot_presence_only.R
#'@date 01.10.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@description plots the results of presence-only metrics

# 1 - install and load packages  ####
#-----------------------------------#

library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(parallel)
library(hexbin)
library(egg)
#library(lmodel2)
library(grid)
library(smatr)


if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

# create a unique name for different runs
nameRun <- paste0("run", 9)

source(paste0("R/",nameRun,"/functions/scale_metric.R"))

# 2 - load data ####
#------------------#

data=readRDS(paste0("data/",nameRun,"/results.RDS"))
data$points <-factor(data$points, levels = c(40, 80, 120, 160, 200, 400))

# List of metrics to scale
metric_names <- c("AUC", "COR", "Spec", "Sens", "Kappa", "SBI_m",
                  "PCC", "TSS", "PRG", "SEDI", "trueCor")

# Apply scaling function to each metric and create new "_scaled" columns
data <- data %>%
  mutate(across(all_of(metric_names),
                .fns = ~ scale_metric(., cur_column()),
                .names = "{.col}_scaled"))
data <- data %>%
  mutate(omissionRate=1-omissionRate,
         SEDI_scaled= 1-SEDI_scaled)

prev=data.frame(species= c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
                prevalence= c(0.35,0.34,0.33,0.29,0.26,0.21,0.15,0.12,0.11,0.05),
                class=c("Broad","Broad","Broad","Middle","Middle","Middle","Narrow","Narrow","Narrow","Narrow"))


#VS1–3 were calibrated to be broadly distributed, VS7–10 narrowly distributed and VS4–6 in between
data=merge(data,prev, by.x="species",by.y="species");rm(prev)


if(!dir.exists(paste0("images/",nameRun,"/resultPlots"))) dir.create(paste0("images/",nameRun,"/resultPlots"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsSummarized"))) dir.create(paste0("images/",nameRun,"/resultPlotsSummarized"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsPrevalence"))) dir.create(paste0("images/",nameRun,"/resultPlotsPrevalence"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsPoint"))) dir.create(paste0("images/",nameRun,"/resultPlotsPoint"), recursive=T)

metric_names <- c("SEDI_scaled", "SBI_m_scaled","omissionRate" )#"SEDI_scaled", "SBI_m_scaled"


labels <-data.frame(metric_names=c("SEDI_scaled", "SBI_m_scaled", "omissionRate"),
                    label = c("SEDI", "'Smoothed boyce index mean'", "'Omission rate'"))


# 3 - plot with reduced major axis regression ####
#------------------------------------------------#

#reduced major axis same as standard major axis (SMA)
# http://stratigrafia.org/8370/lecturenotes/regression.html
# https://stackoverflow.com/questions/49196327/is-there-a-difference-between-reduced-major-axis-regression-and-ranged-major-axi

all_plots=list()
all_results=list()


for (m in metric_names) {
  print(m)
  # for(method_plot in c("PAA", "PA", "PBG")){
  
 # method_plot<- "PA"
  data2=data%>%dplyr::filter(method=="PBG")
  
  # SMA-Regression
  fit <- smatr::sma(paste("trueCor_scaled ~", m),data=data2, method="SMA", slope.test=1, elev.test=0)
  
  # Ergebnisse anzeigen
  #summary(fit)
  
  
  #plot(fit, which="residual" )
  #plot(fit, which="qq")
  
  intercept_rma <- coef(fit)[1][[1]]
  slope_rma <- coef(fit)[2][[1]]
  
  fit_summary<- fit$groupsummary
  
  slope_pValue <- fit_summary$Slope_test_p
  intercept_pValue <- fit_summary$Elev_test_p
  
  # Compute predicted values along the RMA line
  data2$y_rma <- intercept_rma + slope_rma * data2[[m]]
  
  
  # Metrics relative to RMA
  metrics_df <- data2 %>%
    # group_by(method) %>%
    summarise(
      # RMSD: deviation between observed and raw predictor (no regression correction)
      rmsd = sqrt(mean((trueCor_scaled - .data[[m]])^2, na.rm = TRUE)),
      
      # RMSE: residual error after regression correction (with RMA-predicted values)
      #rmse = sqrt(mean((trueCor_scaled - y_rma)^2, na.rm = TRUE)),
      rmse_rma = sqrt(mean((trueCor_scaled - y_rma)^2, na.rm = TRUE)),
      # mae_rma  = mean(abs(trueCor_scaled - y_rma), na.rm = TRUE),
      r_squared_rma = fit_summary$r2,
      #r_squared = rma_fit$rsquare,
      slope_rma=slope_rma,
      intercept_rma=intercept_rma,
      slope_pValue =slope_pValue ,
      intercept_pValue=intercept_pValue,
      # mean_intercept_dev_rma = mean(trueCor_scaled - y_rma, na.rm = TRUE),
      xpos = 0,
      ypos = 0.75,
      .groups = "drop"
    )
  
  # Format significance stars
  slope_sig <- ifelse(slope_pValue < 0.001, "***",
                      ifelse(slope_pValue < 0.01, "**",
                             ifelse(slope_pValue < 0.05, "*", "")))
  intercept_sig <- ifelse(intercept_pValue < 0.001, "***",
                          ifelse(intercept_pValue < 0.01, "**",
                                 ifelse(intercept_pValue < 0.05, "*", "")))
  rsquared_sig <- ifelse(fit_summary$pval < 0.001, "***",
                         ifelse(fit_summary$pval < 0.01, "**",
                                ifelse(fit_summary$pval < 0.05, "*", "")))
  
  label_text <- paste0(
    "RMSE = ", round(metrics_df$rmse_rma, 2),
    "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
    "\nRMSD = ", round(metrics_df$rmsd, 2),
    "\nslope = ", round(slope_rma, 2), slope_sig,
    "\nintercept = ", round(intercept_rma, 2), intercept_sig
  )
  
  # Einmaliger DataFrame für diese Iteration
  annotation_layer_df <- data.frame(
    xpos   = metrics_df$xpos,
    ypos   = metrics_df$ypos,
    label  = label_text
  )
  
  
  p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
    geom_point(size = 0.3, colour="cornflowerblue") +
    geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
    geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
   # facet_wrap(vars(method)) +
    geom_text(data = annotation_layer_df,
              aes(x = xpos, y = ypos, label = label),
              inherit.aes = FALSE, size = 5, hjust = 0) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
    # ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
    ylab("") +
    theme_minimal(base_size = 24)+theme(legend.position="none")
  
 
  
  fit_summary$rmsd <- metrics_df$rmsd
  fit_summary$rmse <- metrics_df$rmse_rma
  fit_summary$intercept <- metrics_df$intercept_rma
  fit_summary$metric <- m
  all_results[[m]]  <- fit_summary
  
  all_plots[[paste0(m)]]  <- p
  
  rm(data2, fit,intercept_rma, slope_rma,p,metrics_df, slope_pValue, intercept_pValue,fit_summary)
  
  
  # }
}

 p1=egg::ggarrange(plots=all_plots,nrow=1,ncol=3, left=textGrob(
  "Pearson correlation between \nsuitability raster and artificial \ndistribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSummarized/presenceOnly_byMethod.png"), dpi = 300, width = 16, height = 5)

