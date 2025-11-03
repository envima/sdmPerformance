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
nameRun <- paste0("run", 13)

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
  
  
  
  # Metrics relative to RMA
  metrics_df <- data2 %>%
    # group_by(method) %>%
    summarise(
      # RMSD: deviation between observed and raw predictor (no regression correction)
      rmse = sqrt(mean((trueCor_scaled - .data[[m]])^2, na.rm = TRUE)),
      xpos = 0,
      ypos = 0.95,
      .groups = "drop"
    )
  
 
  label_text <- paste0(
    "RMSE = ", round(metrics_df$rmse, 2)
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
    #geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
   # facet_wrap(vars(method)) +
    geom_text(data = annotation_layer_df,
              aes(x = xpos, y = ypos, label = label),
              inherit.aes = FALSE, size = 5, hjust = 0) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
    # ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
    ylab("") +
    theme_minimal(base_size = 24)+theme(legend.position="none")
  

  
  all_plots[[paste0(m)]]  <- p
  
  rm(data2, fit,intercept_rma, slope_rma,p,metrics_df, slope_pValue, intercept_pValue,fit_summary)
  
  
  # }
}

 p1=egg::ggarrange(plots=all_plots,nrow=1,ncol=3, left=textGrob(
  "Pearson correlation between \nprobability of occurrence and \nartificial distribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSummarized/presenceOnly_byMethod.png"), dpi = 300, width = 16, height = 5)

