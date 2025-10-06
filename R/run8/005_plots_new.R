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
nameRun <- paste0("run", 8)

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


prev=data.frame(species= c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
                prevalence= c(0.35,0.34,0.33,0.29,0.26,0.21,0.15,0.12,0.11,0.05),
                class=c("Broad","Broad","Broad","Middle","Middle","Middle","Narrow","Narrow","Narrow","Narrow"))


#VS1–3 were calibrated to be broadly distributed, VS7–10 narrowly distributed and VS4–6 in between
data=merge(data,prev, by.x="species",by.y="species");rm(prev)


if(!dir.exists(paste0("images/",nameRun,"/resultPlots"))) dir.create(paste0("images/",nameRun,"/resultPlots"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsSummarized"))) dir.create(paste0("images/",nameRun,"/resultPlotsSummarized"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsPrevalence"))) dir.create(paste0("images/",nameRun,"/resultPlotsPrevalence"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsPoint"))) dir.create(paste0("images/",nameRun,"/resultPlotsPoint"), recursive=T)

metric_names <- c( "AUC_scaled", "COR_scaled", "PRG_scaled", "Spec_scaled", "Sens_scaled", "TSS_scaled", "Kappa_scaled", 
                   "PCC_scaled")#"SEDI_scaled", "SBI_m_scaled"


labels <-data.frame(metric_names=c( "AUC_scaled", "COR_scaled",                              "Spec_scaled", "Sens_scaled", "TSS_scaled", "Kappa_scaled",  "PCC_scaled", "PRG_scaled" ),
                    label = c(        "AUC[ROC]", "'Pearson`s correlation'", "Specificity", "Sensitivity", "TSS",       "'Cohen’s kappa'",    "PCC",     "AUC[PRG]"))



# 3 - plot with reduced major axis regression ####
#------------------------------------------------#

#reduced major axis same as standard major axis (SMA)
# http://stratigrafia.org/8370/lecturenotes/regression.html
# https://stackoverflow.com/questions/49196327/is-there-a-difference-between-reduced-major-axis-regression-and-ranged-major-axi

all_plots=list()
all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    
    data2=data%>%dplyr::filter(method == method_plot)
    
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
      label  = label_text,
      method = method_plot
    )
    
    
    p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
      geom_point(size = 0.3, colour="cornflowerblue") +
      geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
      geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
      facet_wrap(vars(method)) +
      geom_text(data = annotation_layer_df,
                aes(x = xpos, y = ypos, label = label),
                inherit.aes = FALSE, size = 5, hjust = 0) +
      xlim(0, 1) + ylim(-0.01, 1) +
      xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
      # ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
      ylab("") +
      theme_minimal(base_size = 24)+theme(legend.position="none")
    
    assign(paste0("plot_",method_plot), p)
    
    fit_summary$rmsd <- metrics_df$rmsd
    fit_summary$rmse <- metrics_df$rmse_rma
    fit_summary$intercept <- metrics_df$intercept_rma
    fit_summary$metric <- m
    fit_summary$method <- method_plot
    all_results[[paste0(m, "_",method_plot)]]  <- fit_summary
    
    rm(data2, fit,intercept_rma, slope_rma,p,metrics_df, slope_pValue, intercept_pValue,fit_summary)
    
  }
  
  
  
  plots=list(plot_PA, plot_PBG, plot_PAA)
  remove_axis <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  plots[-c(1,4)] <- lapply(plots[-c(1,4)] , function(.p) .p + remove_axis)
  
  p=egg::ggarrange(plots=plots,ncol=3)
  ggsave(p, filename = paste0("images/",nameRun,"/resultPlots/", m, "_byMethod.png"), dpi = 300, width = 16, height = 8)
  all_plots[[paste0(m, "_PA")]]  <- plot_PA
  all_plots[[paste0(m, "_PBG")]] <- plot_PBG
  all_plots[[paste0(m, "_PAA")]] <- plot_PAA
  rm(plot_PA, plot_PBG, plot_PAA,plots,p,remove_axis)
}


#remove_axis <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#all_plots[-c(1,4,7)] <- lapply(all_plots[-c(1,4,7)] , function(.p) .p + remove_axis)
p1=egg::ggarrange(plots=all_plots[1:12],nrow=4,ncol=3, left=textGrob(
  "Pearson correlation between suitability raster and artificial distribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSummarized/plot1_byMethod.png"), dpi = 300, width = 16, height = 20)


p2=egg::ggarrange(plots=all_plots[13:24],nrow=4,ncol=3, left=textGrob(
  "Pearson correlation between suitability raster and artificial distribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsSummarized/plot2_byMethod.png"), dpi = 300, width = 16, height = 20)

rm(p1,p2,labels, slope_sig, rsquared_sig, intercept_sig,method_plot,label_text,all_plots,m)

dataSummary <- do.call(rbind, all_results)
saveRDS(dataSummary, paste0("images/",nameRun,"/resultPlotsSummarized/results.RDS"));rm(dataSummary,all_results)



# 4 - plot with prevalence and reduced major axis regression ####
#---------------------------------------------------------------#


all_plots=list()
all_results=list()

##################################################


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(prevalence_plot in unique(data$class)){
      data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(class == prevalence_plot)
      
      # SMA-Regression
      fit <- smatr::sma(paste("trueCor_scaled ~", m),data=data2, method="SMA", slope.test=1, elev.test=0)
      
      # Ergebnisse anzeigen
      #summary(fit)
      
      
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
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.3, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
        facet_wrap(vars(method,class)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,prevalence_plot), p)
      
           fit_summary$rmsd <- metrics_df$rmsd
      fit_summary$rmse <- metrics_df$rmse_rma
      fit_summary$intercept <- metrics_df$intercept_rma
      fit_summary$metric <- m
      fit_summary$method <- method_plot
      fit_summary$prevalence <- prevalence_plot
      all_results[[paste0(m, "_",method_plot,"_",prevalence_plot)]]  <- fit_summary
      
      rm(data2, fit,intercept_rma,fit_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
    }
  }
  
  
  plots=list(plot_PABroad, plot_PBGBroad, plot_PAABroad,
             plot_PAMiddle, plot_PBGMiddle, plot_PAAMiddle,
             plot_PANarrow, plot_PBGNarrow, plot_PAANarrow)

  p=egg::ggarrange(plots=plots,ncol=3,nrow=3, left=textGrob(
    "Pearson correlation between suitability raster and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p, filename = paste0("images/",nameRun,"/resultPlotsPrevalence/", m, "_byMethodPrevalence.png"), dpi = 300, width = 16, height = 15)
  rm(plot_PABroad, plot_PBGBroad, plot_PAABroad,
     plot_PAMiddle, plot_PBGMiddle, plot_PAAMiddle,
     plot_PANarrow, plot_PBGNarrow, plot_PAANarrow,plots,p)
}

rm(all_plots, data, prev, labels)
dataSummary <- do.call(rbind, all_results)
saveRDS(dataSummary, paste0("images/",nameRun,"/resultPlotsPrevalence/resultsPrevalence.RDS"));rm(dataSummary,all_results)



# 6 - results by number of poiints ####
#-------------------------------------#

all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in c("40", "80", "120", "160", "200", "400")){
      data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(points == point_plot)
      
      # SMA-Regression
      fit <- smatr::sma(paste("trueCor_scaled ~", m),data=data2, method="SMA", slope.test=1, elev.test=0)
      
      # Ergebnisse anzeigen
      #summary(fit)
      
      
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
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.3, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
        facet_wrap(vars(method,points)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,point_plot), p)
      
      fit_summary$rmsd <- metrics_df$rmsd
      fit_summary$rmse <- metrics_df$rmse_rma
      fit_summary$intercept <- metrics_df$intercept_rma
      fit_summary$metric <- m
      fit_summary$method <- method_plot
      fit_summary$point <- point_plot
      all_results[[paste0(m, "_",method_plot,"_",point_plot)]]  <- fit_summary
      
      rm(data2, fit,intercept_rma,fit_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
    }
  }
  
  
  plots=list(plot_PA40,plot_PBG40,plot_PAA40,
             plot_PA80,plot_PBG80,plot_PAA80,
             plot_PA120,plot_PBG120,plot_PAA120,
             plot_PA160,plot_PBG160,plot_PAA160,
             plot_PA200,plot_PBG200,plot_PAA200,
             plot_PA400,plot_PBG400,plot_PAA400)
  rm(plot_PA40,plot_PBG40,plot_PAA40,plot_PA80,plot_PBG80,plot_PAA80,plot_PA120,plot_PBG120,plot_PAA120,plot_PA160,plot_PBG160,plot_PAA160,plot_PA200,plot_PBG200,plot_PAA200,plot_PA400,plot_PBG400,plot_PAA400)
  
  p1=egg::ggarrange(plots=plots[1:12],ncol=3,nrow=4, left=textGrob(
    "Pearson correlation between suitability raster /nand artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsPoint/", m, "_byMethodPoint1.png"), dpi = 300, width = 16, height = 20)
  
  p2=egg::ggarrange(plots=plots[13:18],ncol=3,nrow=2, left=textGrob(
    "Pearson correlation between suitability \nraster and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsPoint/", m, "_byMethodPoint2.png"), dpi = 300, width = 16, height = 10)
  
  
  
  rm(plots,p1,p2)
}

rm(all_plots, data, prev, labels,annotation_layer_df)
dataSummary <- do.call(rbind, all_results)
saveRDS(dataSummary, paste0("images/",nameRun,"/resultPlotsPoint/resultsPoint.RDS"))
write.csv(dataSummary, paste0("images/",nameRun,"/resultPlotsPoint/resultsPoint.csv"));rm(dataSummary,all_results)
