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
#library(smatr)
#library(SimplyAgree)
#library("mcr")

if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

# create a unique name for different runs
nameRun <- paste0("run", 14)

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

if(!dir.exists(paste0("images/",nameRun,"/resultPlotsModel"))) dir.create(paste0("images/",nameRun,"/resultPlotsModel"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsSize"))) dir.create(paste0("images/",nameRun,"/resultPlotsSize"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsSizeModel"))) dir.create(paste0("images/",nameRun,"/resultPlotsSizeModel"), recursive=T)
if(!dir.exists(paste0("images/",nameRun,"/resultPlotsSizeModelPoints"))) dir.create(paste0("images/",nameRun,"/resultPlotsSizeModelPoints"), recursive=T)




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
    #PBequi <- mcreg(data2[["trueCor_scaled"]],data2[[m]], 
    #                method.reg = "PBequi", 
    #                slope.measure="radian", 
    #                method.ci="analytical",
    #                methodlarge=T)
    #PBequi@para
    
    
    #mcr::calcCUSUM(PBequi)
    
    #summary(PBequi)
    #intercept <- coef(PBequi)[1][[1]]
    #slope <- coef(PBequi)[2][[1]]
    
    #PBequi_summary<-  summary(PBequi)
    
    # Compute predicted values along the RMA line
    # data2$y_PBequi <- intercept + slope * data2[[m]]
    #  
    # Metrics relative to RMA
    metrics_df <- data2 %>%
      # group_by(method) %>%
      summarise(
        # RMSD: deviation between observed and raw predictor (no regression correction)
        rmse = sqrt(mean((trueCor_scaled - .data[[m]])^2, na.rm = TRUE)),
        # RMSE: residual error after regression correction (with RMA-predicted values)
        #rmse2 = sqrt(mean((trueCor_scaled - y_PBequi)^2, na.rm = TRUE)),
        #rmse_rma = sqrt(mean((trueCor_scaled - y_rma)^2, na.rm = TRUE)),
        #mae_rma  = mean(abs(trueCor_scaled - y_rma), na.rm = TRUE),
        #r_squared = rma_PBequi$rsquare,
        # slope=slope,
        #intercept=intercept,
        #mean_intercept_dev_rma = mean(trueCor_scaled - y_rma, na.rm = TRUE),
        xpos = 0,
        ypos = 0.95,
        .groups = "drop"
      )
    
    ## Format significance stars
    #slope_sig <- ifelse(slope_pValue < 0.001, "***",
    #                    ifelse(slope_pValue < 0.01, "**",
    #                           ifelse(slope_pValue < 0.05, "*", "")))
    #intercept_sig <- ifelse(intercept_pValue < 0.001, "***",
    #                        ifelse(intercept_pValue < 0.01, "**",
    #                               ifelse(intercept_pValue < 0.05, "*", "")))
    #rsquared_sig <- ifelse(PBequi_summary$pval < 0.001, "***",
    #                       ifelse(PBequi_summary$pval < 0.01, "**",
    #                              ifelse(PBequi_summary$pval < 0.05, "*", "")))
    
    label_text <- paste0(
      # "RMSE2 = ", round(metrics_df$rmse2, 2),
      #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
      "\nRMSE = ", round(metrics_df$rmse, 2)
      # "\nslope = ", round(slope, 2),
      # "\nintercept = ", round(intercept, 2)
    )
    
    # Einmaliger DataFrame für diese Iteration
    annotation_layer_df <- data.frame(
      xpos   = metrics_df$xpos,
      ypos   = metrics_df$ypos,
      label  = label_text,
      method = method_plot
    )
    
    
    p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
      geom_point(size = 0.1, colour="cornflowerblue") +
      geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
      # geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
      
      facet_wrap(vars(method)) +
      geom_text(data = annotation_layer_df,
                aes(x = xpos, y = ypos, label = label),
                inherit.aes = FALSE, size = 5, hjust = 0) +
      xlim(0, 1) + ylim(-0.01, 1) +
      xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
      # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
      ylab("") +
      theme_minimal(base_size = 24)+theme(legend.position="none")
    
    assign(paste0("plot_",method_plot), p)
    
    # PBequi_summary$rmsd <- metrics_df$rmsd
    #  PBequi_summary$rmse <- metrics_df$rmse_rma
    #  PBequi_summary$intercept <- metrics_df$intercept_rma
    #  PBequi_summary$metric <- m
    #  PBequi_summary$method <- method_plot
    #  all_results[[paste0(m, "_",method_plot)]]  <- PBequi_summary
    
    rm(data2, PBequi,intercept, slope,p,metrics_df)
    
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
  "Pearson correlation between probability of occurrence and artificial distribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSummarized/plot1_byMethod.png"), dpi = 300, width = 16, height = 20)


p2=egg::ggarrange(plots=all_plots[13:24],nrow=4,ncol=3, left=textGrob(
  "Pearson correlation between probability of occurrence and artificial distribution map",
  gp = gpar(fontsize = 24),
  rot = 90))
ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsSummarized/plot2_byMethod.png"), dpi = 300, width = 16, height = 20)

rm(p1,p2, method_plot,label_text,all_plots,m)

#dataSummary <- do.call(rbind, all_results)
#saveRDS(dataSummary, paste0("images/",nameRun,"/resultPlotsSummarized/results.RDS"));rm(dataSummary,all_results)



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
        "\nRMSE = ", round(metrics_df$rmse, 2)
      )
      
      # Einmaliger DataFrame für diese Iteration
      annotation_layer_df <- data.frame(
        xpos   = metrics_df$xpos,
        ypos   = metrics_df$ypos,
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.1, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        #geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
        facet_wrap(vars(method,class)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,prevalence_plot), p)
      
      
      rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
    }
  }
  
  
  plots=list(plot_PABroad, plot_PBGBroad, plot_PAABroad,
             plot_PAMiddle, plot_PBGMiddle, plot_PAAMiddle,
             plot_PANarrow, plot_PBGNarrow, plot_PAANarrow)
  
  p=egg::ggarrange(plots=plots,ncol=3,nrow=3, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p, filename = paste0("images/",nameRun,"/resultPlotsPrevalence/", m, "_byMethodPrevalence.png"), dpi = 300, width = 16, height = 15)
  rm(plot_PABroad, plot_PBGBroad, plot_PAABroad,
     plot_PAMiddle, plot_PBGMiddle, plot_PAAMiddle,
     plot_PANarrow, plot_PBGNarrow, plot_PAANarrow,plots,p)
}




# 6 - results by number of poiints ####
#-------------------------------------#

all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in c("40", "80", "120", "160", "200", "400")){
      data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(points == point_plot)
      
      
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
        #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
        "\nRMSE = ", round(metrics_df$rmse, 2)
      )
      
      
      # Einmaliger DataFrame für diese Iteration
      annotation_layer_df <- data.frame(
        xpos   = metrics_df$xpos,
        ypos   = metrics_df$ypos,
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.1, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        #  geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
        facet_wrap(vars(method,points)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,point_plot), p)
      
      
      
      rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
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
    "Pearson correlation between probability of \noccurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsPoint/", m, "_byMethodPoint1.png"), dpi = 300, width = 16, height = 20)
  
  p2=egg::ggarrange(plots=plots[13:18],ncol=3,nrow=2, left=textGrob(
    "Pearson correlation between probability of \noccurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsPoint/", m, "_byMethodPoint2.png"), dpi = 300, width = 16, height = 10)
  
  
  
  rm(plots,p1,p2)
}

rm(all_plots, data,  labels,annotation_layer_df)


# 7 - results by model ####
#-------------------------#

all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in unique(data$model)){
      data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(model == point_plot)
      
      
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
        #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
        "\nRMSE = ", round(metrics_df$rmse, 2)
      )
      
      
      # Einmaliger DataFrame für diese Iteration
      annotation_layer_df <- data.frame(
        xpos   = metrics_df$xpos,
        ypos   = metrics_df$ypos,
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.1, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        #  geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
        facet_wrap(vars(method,model)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,"_",point_plot), p)
      
      
      
      rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
    }
  }
  
  
  plots=list(plot_PA_BRT,plot_PBG_BRT,plot_PAA_BRT,
             plot_PA_GAM,plot_PBG_GAM,plot_PAA_GAM,
             plot_PA_Maxent,plot_PBG_Maxent,plot_PAA_Maxent,
             plot_PA_Lasso,plot_PBG_Lasso,plot_PAA_Lasso,
             plot_PA_RF,plot_PBG_RF,plot_PAA_RF)
  
  rm(plot_PA_BRT,plot_PBG_BRT,plot_PAA_BRT,
     plot_PA_GAM,plot_PBG_GAM,plot_PAA_GAM,
     plot_PA_Maxent,plot_PBG_Maxent,plot_PAA_Maxent,
     plot_PA_Lasso,plot_PBG_Lasso,plot_PAA_Lasso,
     plot_PA_RF,plot_PBG_RF,plot_PAA_RF)
  
  p1=egg::ggarrange(plots=plots[1:15],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsModel/", m, "_byMethodModel.png"), dpi = 300, width = 16, height = 25)
  
  rm(plots,p1)
}

rm(all_plots, data,  labels,annotation_layer_df)


# 7 - results by size ####
#------------------------#

all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in unique(data$size)){
      data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(size == point_plot)
      
      
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
        #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
        "\nRMSE = ", round(metrics_df$rmse, 2)
      )
      
      
      # Einmaliger DataFrame für diese Iteration
      annotation_layer_df <- data.frame(
        xpos   = metrics_df$xpos,
        ypos   = metrics_df$ypos,
        label  = label_text,
        method = method_plot
      )
      
      
      p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
        geom_point(size = 0.1, colour="cornflowerblue") +
        geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
        #  geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
        facet_wrap(vars(method,size)) +
        geom_text(data = annotation_layer_df,
                  aes(x = xpos, y = ypos, label = label),
                  inherit.aes = FALSE, size = 5, hjust = 0) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlim(0, 1) + ylim(-0.01, 1) +
        xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
        # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
        ylab("") +
        theme_minimal(base_size = 24)+theme(legend.position="none")
      
      assign(paste0("plot_",method_plot,"_",point_plot), p)
      
      
      
      rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      
    }
  }
  
  
  plots=list(plot_PA_block1,plot_PBG_block1,plot_PAA_block1,
             plot_PA_block2,plot_PBG_block2,plot_PAA_block2,
             plot_PA_clusters,plot_PBG_clusters,plot_PAA_clusters,
             plot_PA_KNNDM,plot_PBG_KNNDM,plot_PAA_KNNDM,
             plot_PA_random,plot_PBG_random,plot_PAA_random)
  
  rm(plot_PA_block1,plot_PBG_block1,plot_PAA_block1,
     plot_PA_block2,plot_PBG_block2,plot_PAA_block2,
     plot_PA_clusters,plot_PBG_clusters,plot_PAA_clusters,
     plot_PA_KNNDM,plot_PBG_KNNDM,plot_PAA_KNNDM,
     plot_PA_random,plot_PBG_random,plot_PAA_random)
  
  p1=egg::ggarrange(plots=plots[1:15],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSize/", m, "_byMethodSize.png"), dpi = 300, width = 16, height = 25)
  
  rm(plots,p1)
}


# 8 - results by size and model ####
#----------------------------------#

all_results=list()


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in unique(data$size)){
      for(model_plot in unique(data$model)){
        data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(size == point_plot)%>%dplyr::filter(model == model_plot)
        
        
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
          #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
          "\nRMSE = ", round(metrics_df$rmse, 2)
        )
        
        
        # Einmaliger DataFrame für diese Iteration
        annotation_layer_df <- data.frame(
          xpos   = metrics_df$xpos,
          ypos   = metrics_df$ypos,
          label  = label_text,
          method = method_plot
        )
        
        
        p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
          geom_point(size = 0.1, colour="cornflowerblue") +
          geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
          #  geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
          facet_wrap(vars(method,size,model)) +
          geom_text(data = annotation_layer_df,
                    aes(x = xpos, y = ypos, label = label),
                    inherit.aes = FALSE, size = 5, hjust = 0) +
          xlim(0, 1) + ylim(-0.01, 1) +
          xlim(0, 1) + ylim(-0.01, 1) +
          xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
          # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
          ylab("") +
          theme_minimal(base_size = 24)+theme(legend.position="none")
        
        assign(paste0("plot_",method_plot,"_",point_plot,"_",model_plot), p)
        
        
        
        rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
      }
    }
    
    
    
    
  }
  
  
  plots=list(plot_PA_block1_BRT,plot_PBG_block1_BRT,plot_PAA_block1_BRT,
             plot_PA_block2_BRT,plot_PBG_block2_BRT,plot_PAA_block2_BRT,
             plot_PA_clusters_BRT,plot_PBG_clusters_BRT,plot_PAA_clusters_BRT,
             plot_PA_KNNDM_BRT,plot_PBG_KNNDM_BRT,plot_PAA_KNNDM_BRT,
             plot_PA_random_BRT,plot_PBG_random_BRT,plot_PAA_random_BRT,
             plot_PA_block1_GAM,plot_PBG_block1_GAM,plot_PAA_block1_GAM,
             plot_PA_block2_GAM,plot_PBG_block2_GAM,plot_PAA_block2_GAM,
             plot_PA_clusters_GAM,plot_PBG_clusters_GAM,plot_PAA_clusters_GAM,
             plot_PA_KNNDM_GAM,plot_PBG_KNNDM_GAM,plot_PAA_KNNDM_GAM,
             plot_PA_random_GAM,plot_PBG_random_GAM,plot_PAA_random_GAM,
             plot_PA_block1_Lasso,plot_PBG_block1_Lasso,plot_PAA_block1_Lasso,
             plot_PA_block2_Lasso,plot_PBG_block2_Lasso,plot_PAA_block2_Lasso,
             plot_PA_clusters_Lasso,plot_PBG_clusters_Lasso,plot_PAA_clusters_Lasso,
             plot_PA_KNNDM_Lasso,plot_PBG_KNNDM_Lasso,plot_PAA_KNNDM_Lasso,
             plot_PA_random_Lasso,plot_PBG_random_Lasso,plot_PAA_random_Lasso,
             plot_PA_block1_Maxent,plot_PBG_block1_Maxent,plot_PAA_block1_Maxent,
             plot_PA_block2_Maxent,plot_PBG_block2_Maxent,plot_PAA_block2_Maxent,
             plot_PA_clusters_Maxent,plot_PBG_clusters_Maxent,plot_PAA_clusters_Maxent,
             plot_PA_KNNDM_Maxent,plot_PBG_KNNDM_Maxent,plot_PAA_KNNDM_Maxent,
             plot_PA_random_Maxent,plot_PBG_random_Maxent,plot_PAA_random_Maxent,
             plot_PA_block1_RF,plot_PBG_block1_RF,plot_PAA_block1_RF,
             plot_PA_block2_RF,plot_PBG_block2_RF,plot_PAA_block2_RF,
             plot_PA_clusters_RF,plot_PBG_clusters_RF,plot_PAA_clusters_RF,
             plot_PA_KNNDM_RF,plot_PBG_KNNDM_RF,plot_PAA_KNNDM_RF,
             plot_PA_random_RF,plot_PBG_random_RF,plot_PAA_random_RF)
  
  rm(plot_PA_block1_BRT,plot_PBG_block1_BRT,plot_PAA_block1_BRT,
     plot_PA_block2_BRT,plot_PBG_block2_BRT,plot_PAA_block2_BRT,
     plot_PA_clusters_BRT,plot_PBG_clusters_BRT,plot_PAA_clusters_BRT,
     plot_PA_KNNDM_BRT,plot_PBG_KNNDM_BRT,plot_PAA_KNNDM_BRT,  
     plot_PA_random_BRT,plot_PBG_random_BRT,plot_PAA_random_BRT,
     plot_PA_block1_GAM,plot_PBG_block1_GAM,plot_PAA_block1_GAM,
     plot_PA_block2_GAM,plot_PBG_block2_GAM,plot_PAA_block2_GAM,
     plot_PA_clusters_GAM,plot_PBG_clusters_GAM,plot_PAA_clusters_GAM,
     plot_PA_KNNDM_GAM,plot_PBG_KNNDM_GAM,plot_PAA_KNNDM_GAM,
     plot_PA_random_GAM,plot_PBG_random_GAM,plot_PAA_random_GAM,
     plot_PA_block1_Lasso,plot_PBG_block1_Lasso,plot_PAA_block1_Lasso,
     plot_PA_block2_Lasso,plot_PBG_block2_Lasso,plot_PAA_block2_Lasso,
     plot_PA_clusters_Lasso,plot_PBG_clusters_Lasso,plot_PAA_clusters_Lasso,
     plot_PA_KNNDM_Lasso,plot_PBG_KNNDM_Lasso,plot_PAA_KNNDM_Lasso,
     plot_PA_random_Lasso,plot_PBG_random_Lasso,plot_PAA_random_Lasso,
     plot_PA_block1_Maxent,plot_PBG_block1_Maxent,plot_PAA_block1_Maxent,
     plot_PA_block2_Maxent,plot_PBG_block2_Maxent,plot_PAA_block2_Maxent,
     plot_PA_clusters_Maxent,plot_PBG_clusters_Maxent,plot_PAA_clusters_Maxent,
     plot_PA_KNNDM_Maxent,plot_PBG_KNNDM_Maxent,plot_PAA_KNNDM_Maxent,
     plot_PA_random_Maxent,plot_PBG_random_Maxent,plot_PAA_random_Maxent,
     plot_PA_block1_RF,plot_PBG_block1_RF,plot_PAA_block1_RF,
     plot_PA_block2_RF,plot_PBG_block2_RF,plot_PAA_block2_RF,
     plot_PA_clusters_RF,plot_PBG_clusters_RF,plot_PAA_clusters_RF,
     plot_PA_KNNDM_RF,plot_PBG_KNNDM_RF,plot_PAA_KNNDM_RF,
     plot_PA_random_RF,plot_PBG_random_RF,plot_PAA_random_RF)
  
  p1=egg::ggarrange(plots=plots[1:15],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSizeModel/", m, "_BRT_byMethodSizeModel.png"), dpi = 300, width = 16, height = 25)
  
  p2=egg::ggarrange(plots=plots[16:30],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsSizeModel/", m, "_GAM_byMethodSizeModel.png"), dpi = 300, width = 16, height = 25)
  
  p3=egg::ggarrange(plots=plots[31:45],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p3, filename = paste0("images/",nameRun,"/resultPlotsSizeModel/", m, "_Lasso_byMethodSizeModel.png"), dpi = 300, width = 16, height = 25)
  
  p4=egg::ggarrange(plots=plots[46:60],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p4, filename = paste0("images/",nameRun,"/resultPlotsSizeModel/", m, "_Maxent_byMethodSizeModel.png"), dpi = 300, width = 16, height = 25)
  
  
  p5=egg::ggarrange(plots=plots[61:75],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p5, filename = paste0("images/",nameRun,"/resultPlotsSizeModel/", m, "_RF_byMethodSizeModel.png"), dpi = 300, width = 16, height = 25)
  
  
  
  rm(plots,p1,p2,p3,p4,p5)
}

# 9 - results by size and model and points ####
#---------------------------------------------#

all_results=list()

metric_names <- c( "AUC_scaled", "COR_scaled",  "TSS_scaled" )
                    #"SEDI_scaled", "SBI_m_scaled","Spec_scaled", "Sens_scaled","Kappa_scaled", "PCC_scaled","PRG_scaled", 


for (m in metric_names) {
  print(m)
  for(method_plot in c("PAA", "PA", "PBG")){
    for(point_plot in unique(data$points)){
      for(model_plot in unique(data$model)){
        for(size_plot in unique(data$size)){
          data2=data%>%dplyr::filter(method == method_plot)%>%dplyr::filter(points == point_plot)%>%dplyr::filter(model == model_plot)%>%dplyr::filter(size == size_plot)
          
          
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
            #   "\nR² = ", round(metrics_df$r_squared_rma, 2),rsquared_sig,
            "\nRMSE = ", round(metrics_df$rmse, 2)
          )
          
          
          # Einmaliger DataFrame für diese Iteration
          annotation_layer_df <- data.frame(
            xpos   = metrics_df$xpos,
            ypos   = metrics_df$ypos,
            label  = label_text,
            method = method_plot
          )
          
          
          p<-  ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
            geom_point(size = 0.6, colour="cornflowerblue") +
            geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
            #  geom_abline(slope = slope, intercept = intercept, linewidth=1) + # RMA line
            facet_wrap(vars(method,size,model, points)) +
            geom_text(data = annotation_layer_df,
                      aes(x = xpos, y = ypos, label = label),
                      inherit.aes = FALSE, size = 5, hjust = 0) +
            xlim(0, 1) + ylim(-0.01, 1) +
            xlim(0, 1) + ylim(-0.01, 1) +
            xlab(parse(text =  labels[labels$metric_names == m,]$label)) +
            # ylab("Pearson correlation between probability of \noccurrence and artificial distribution map") +
            ylab("") +
            theme_minimal(base_size = 24)+theme(legend.position="none")
          
          assign(paste0("plot_",method_plot,"_",point_plot,"_", size_plot, "_", model_plot), p)
          
          
          
          rm(data2, PBequi,intercept_rma,PBequi_summary, slope_rma,p,metrics_df, slope_pValue, intercept_pValue, intercept_sig, label_text, rsquared_sig, slope_sig)
        }
      }
    } 
    
    
    
  }
  
  
  plots=list(plot_PA_40_block1_BRT,      plot_PBG_40_block1_BRT,              plot_PAA_40_block1_BRT,
             plot_PA_40_block2_BRT,      plot_PBG_40_block2_BRT,              plot_PAA_40_block2_BRT,
             plot_PA_40_clusters_BRT,    plot_PBG_40_clusters_BRT,            plot_PAA_40_clusters_BRT,
             plot_PA_40_KNNDM_BRT,       plot_PBG_40_KNNDM_BRT,               plot_PAA_40_KNNDM_BRT,
             plot_PA_40_random_BRT,      plot_PBG_40_random_BRT,              plot_PAA_40_random_BRT,
             plot_PA_40_block1_GAM,      plot_PBG_40_block1_GAM,              plot_PAA_40_block1_GAM,
             plot_PA_40_block2_GAM,      plot_PBG_40_block2_GAM,              plot_PAA_40_block2_GAM,
             plot_PA_40_clusters_GAM,    plot_PBG_40_clusters_GAM,            plot_PAA_40_clusters_GAM,
             plot_PA_40_KNNDM_GAM,       plot_PBG_40_KNNDM_GAM,               plot_PAA_40_KNNDM_GAM,
             plot_PA_40_random_GAM,      plot_PBG_40_random_GAM,              plot_PAA_40_random_GAM,
             plot_PA_40_block1_Lasso,    plot_PBG_40_block1_Lasso,            plot_PAA_40_block1_Lasso,
             plot_PA_40_block2_Lasso,    plot_PBG_40_block2_Lasso,            plot_PAA_40_block2_Lasso,
             plot_PA_40_clusters_Lasso,  plot_PBG_40_clusters_Lasso,          plot_PAA_40_clusters_Lasso,
             plot_PA_40_KNNDM_Lasso,     plot_PBG_40_KNNDM_Lasso,             plot_PAA_40_KNNDM_Lasso,
             plot_PA_40_random_Lasso,    plot_PBG_40_random_Lasso,            plot_PAA_40_random_Lasso,
             plot_PA_40_block1_Maxent,   plot_PBG_40_block1_Maxent,           plot_PAA_40_block1_Maxent,
             plot_PA_40_block2_Maxent,   plot_PBG_40_block2_Maxent,           plot_PAA_40_block2_Maxent,
             plot_PA_40_clusters_Maxent, plot_PBG_40_clusters_Maxent,         plot_PAA_40_clusters_Maxent,
             plot_PA_40_KNNDM_Maxent,    plot_PBG_40_KNNDM_Maxent,            plot_PAA_40_KNNDM_Maxent,
             plot_PA_40_random_Maxent,   plot_PBG_40_random_Maxent,           plot_PAA_40_random_Maxent,
             plot_PA_40_block1_RF,       plot_PBG_40_block1_RF,               plot_PAA_40_block1_RF,
             plot_PA_40_block2_RF,       plot_PBG_40_block2_RF,               plot_PAA_40_block2_RF,
             plot_PA_40_clusters_RF,     plot_PBG_40_clusters_RF,             plot_PAA_40_clusters_RF,
             plot_PA_40_KNNDM_RF,        plot_PBG_40_KNNDM_RF,                plot_PAA_40_KNNDM_RF,
             plot_PA_40_random_RF,       plot_PBG_40_random_RF,               plot_PAA_40_random_RF,
             plot_PA_80_block1_BRT,      plot_PBG_80_block1_BRT,              plot_PAA_80_block1_BRT,
             plot_PA_80_block2_BRT,      plot_PBG_80_block2_BRT,              plot_PAA_80_block2_BRT,
             plot_PA_80_clusters_BRT,    plot_PBG_80_clusters_BRT,            plot_PAA_80_clusters_BRT,
             plot_PA_80_KNNDM_BRT,       plot_PBG_80_KNNDM_BRT,               plot_PAA_80_KNNDM_BRT,
             plot_PA_80_random_BRT,      plot_PBG_80_random_BRT,              plot_PAA_80_random_BRT,
             plot_PA_80_block1_GAM,      plot_PBG_80_block1_GAM,              plot_PAA_80_block1_GAM,
             plot_PA_80_block2_GAM,      plot_PBG_80_block2_GAM,              plot_PAA_80_block2_GAM,
             plot_PA_80_clusters_GAM,    plot_PBG_80_clusters_GAM,            plot_PAA_80_clusters_GAM,
             plot_PA_80_KNNDM_GAM,       plot_PBG_80_KNNDM_GAM,               plot_PAA_80_KNNDM_GAM,
             plot_PA_80_random_GAM,      plot_PBG_80_random_GAM,              plot_PAA_80_random_GAM,
             plot_PA_80_block1_Lasso,    plot_PBG_80_block1_Lasso,            plot_PAA_80_block1_Lasso,
             plot_PA_80_block2_Lasso,    plot_PBG_80_block2_Lasso,            plot_PAA_80_block2_Lasso,
             plot_PA_80_clusters_Lasso,  plot_PBG_80_clusters_Lasso,          plot_PAA_80_clusters_Lasso,
             plot_PA_80_KNNDM_Lasso,     plot_PBG_80_KNNDM_Lasso,             plot_PAA_80_KNNDM_Lasso,
             plot_PA_80_random_Lasso,    plot_PBG_80_random_Lasso,            plot_PAA_80_random_Lasso,
             plot_PA_80_block1_Maxent,   plot_PBG_80_block1_Maxent,           plot_PAA_80_block1_Maxent,
             plot_PA_80_block2_Maxent,   plot_PBG_80_block2_Maxent,           plot_PAA_80_block2_Maxent,
             plot_PA_80_clusters_Maxent, plot_PBG_80_clusters_Maxent,         plot_PAA_80_clusters_Maxent,
             plot_PA_80_KNNDM_Maxent,    plot_PBG_80_KNNDM_Maxent,            plot_PAA_80_KNNDM_Maxent,
             plot_PA_80_random_Maxent,   plot_PBG_80_random_Maxent,           plot_PAA_80_random_Maxent,
             plot_PA_80_block1_RF,       plot_PBG_80_block1_RF,               plot_PAA_80_block1_RF,
             plot_PA_80_block2_RF,       plot_PBG_80_block2_RF,               plot_PAA_80_block2_RF,
             plot_PA_80_clusters_RF,     plot_PBG_80_clusters_RF,             plot_PAA_80_clusters_RF,
             plot_PA_80_KNNDM_RF,        plot_PBG_80_KNNDM_RF,                plot_PAA_80_KNNDM_RF,
             plot_PA_80_random_RF,       plot_PBG_80_random_RF,               plot_PAA_80_random_RF,
             plot_PA_120_block1_BRT,     plot_PBG_120_block1_BRT,             plot_PAA_120_block1_BRT,
             plot_PA_120_block2_BRT,     plot_PBG_120_block2_BRT,             plot_PAA_120_block2_BRT,
             plot_PA_120_clusters_BRT,   plot_PBG_120_clusters_BRT,           plot_PAA_120_clusters_BRT,
             plot_PA_120_KNNDM_BRT,      plot_PBG_120_KNNDM_BRT,              plot_PAA_120_KNNDM_BRT,
             plot_PA_120_random_BRT,     plot_PBG_120_random_BRT,             plot_PAA_120_random_BRT,
             plot_PA_120_block1_GAM,     plot_PBG_120_block1_GAM,             plot_PAA_120_block1_GAM,
             plot_PA_120_block2_GAM,     plot_PBG_120_block2_GAM,             plot_PAA_120_block2_GAM,
             plot_PA_120_clusters_GAM,   plot_PBG_120_clusters_GAM,           plot_PAA_120_clusters_GAM,
             plot_PA_120_KNNDM_GAM,      plot_PBG_120_KNNDM_GAM,              plot_PAA_120_KNNDM_GAM,
             plot_PA_120_random_GAM,     plot_PBG_120_random_GAM,             plot_PAA_120_random_GAM,
             plot_PA_120_block1_Lasso,   plot_PBG_120_block1_Lasso,           plot_PAA_120_block1_Lasso,
             plot_PA_120_block2_Lasso,   plot_PBG_120_block2_Lasso,           plot_PAA_120_block2_Lasso,
             plot_PA_120_clusters_Lasso, plot_PBG_120_clusters_Lasso,         plot_PAA_120_clusters_Lasso,
             plot_PA_120_KNNDM_Lasso,    plot_PBG_120_KNNDM_Lasso,            plot_PAA_120_KNNDM_Lasso,
             plot_PA_120_random_Lasso,   plot_PBG_120_random_Lasso,           plot_PAA_120_random_Lasso,
             plot_PA_120_block1_Maxent,  plot_PBG_120_block1_Maxent,          plot_PAA_120_block1_Maxent,
             plot_PA_120_block2_Maxent,  plot_PBG_120_block2_Maxent,          plot_PAA_120_block2_Maxent,
             plot_PA_120_clusters_Maxent,plot_PBG_120_clusters_Maxent,        plot_PAA_120_clusters_Maxent,
             plot_PA_120_KNNDM_Maxent,   plot_PBG_120_KNNDM_Maxent,           plot_PAA_120_KNNDM_Maxent,
             plot_PA_120_random_Maxent,  plot_PBG_120_random_Maxent,          plot_PAA_120_random_Maxent,
             plot_PA_120_block1_RF,      plot_PBG_120_block1_RF,              plot_PAA_120_block1_RF,
             plot_PA_120_block2_RF,      plot_PBG_120_block2_RF,              plot_PAA_120_block2_RF,
             plot_PA_120_clusters_RF,    plot_PBG_120_clusters_RF,            plot_PAA_120_clusters_RF,
             plot_PA_120_KNNDM_RF,       plot_PBG_120_KNNDM_RF,               plot_PAA_120_KNNDM_RF,
             plot_PA_120_random_RF,      plot_PBG_120_random_RF,              plot_PAA_120_random_RF,
             plot_PA_160_block1_BRT,     plot_PBG_160_block1_BRT,             plot_PAA_160_block1_BRT,
             plot_PA_160_block2_BRT,     plot_PBG_160_block2_BRT,             plot_PAA_160_block2_BRT,
             plot_PA_160_clusters_BRT,   plot_PBG_160_clusters_BRT,           plot_PAA_160_clusters_BRT,
             plot_PA_160_KNNDM_BRT,      plot_PBG_160_KNNDM_BRT,              plot_PAA_160_KNNDM_BRT,
             plot_PA_160_random_BRT,     plot_PBG_160_random_BRT,             plot_PAA_160_random_BRT,
             plot_PA_160_block1_GAM,     plot_PBG_160_block1_GAM,             plot_PAA_160_block1_GAM,
             plot_PA_160_block2_GAM,     plot_PBG_160_block2_GAM,             plot_PAA_160_block2_GAM,
             plot_PA_160_clusters_GAM,   plot_PBG_160_clusters_GAM,           plot_PAA_160_clusters_GAM,
             plot_PA_160_KNNDM_GAM,      plot_PBG_160_KNNDM_GAM,              plot_PAA_160_KNNDM_GAM,
             plot_PA_160_random_GAM,     plot_PBG_160_random_GAM,             plot_PAA_160_random_GAM,
             plot_PA_160_block1_Lasso,   plot_PBG_160_block1_Lasso,           plot_PAA_160_block1_Lasso,
             plot_PA_160_block2_Lasso,   plot_PBG_160_block2_Lasso,           plot_PAA_160_block2_Lasso,
             plot_PA_160_clusters_Lasso, plot_PBG_160_clusters_Lasso,         plot_PAA_160_clusters_Lasso,
             plot_PA_160_KNNDM_Lasso,    plot_PBG_160_KNNDM_Lasso,            plot_PAA_160_KNNDM_Lasso,
             plot_PA_160_random_Lasso,   plot_PBG_160_random_Lasso,           plot_PAA_160_random_Lasso,
             plot_PA_160_block1_Maxent,  plot_PBG_160_block1_Maxent,          plot_PAA_160_block1_Maxent,
             plot_PA_160_block2_Maxent,  plot_PBG_160_block2_Maxent,          plot_PAA_160_block2_Maxent,
             plot_PA_160_clusters_Maxent,plot_PBG_160_clusters_Maxent,        plot_PAA_160_clusters_Maxent,
             plot_PA_160_KNNDM_Maxent,   plot_PBG_160_KNNDM_Maxent,           plot_PAA_160_KNNDM_Maxent,
             plot_PA_160_random_Maxent,  plot_PBG_160_random_Maxent,          plot_PAA_160_random_Maxent,
             plot_PA_160_block1_RF,      plot_PBG_160_block1_RF,              plot_PAA_160_block1_RF,
             plot_PA_160_block2_RF,      plot_PBG_160_block2_RF,              plot_PAA_160_block2_RF,
             plot_PA_160_clusters_RF,    plot_PBG_160_clusters_RF,            plot_PAA_160_clusters_RF,
             plot_PA_160_KNNDM_RF,       plot_PBG_160_KNNDM_RF,               plot_PAA_160_KNNDM_RF,
             plot_PA_160_random_RF,      plot_PBG_160_random_RF,              plot_PAA_160_random_RF,
             plot_PA_200_block1_BRT,     plot_PBG_200_block1_BRT,             plot_PAA_200_block1_BRT,
             plot_PA_200_block2_BRT,     plot_PBG_200_block2_BRT,             plot_PAA_200_block2_BRT,
             plot_PA_200_clusters_BRT,   plot_PBG_200_clusters_BRT,           plot_PAA_200_clusters_BRT,
             plot_PA_200_KNNDM_BRT,      plot_PBG_200_KNNDM_BRT,              plot_PAA_200_KNNDM_BRT,
             plot_PA_200_random_BRT,     plot_PBG_200_random_BRT,             plot_PAA_200_random_BRT,
             plot_PA_200_block1_GAM,     plot_PBG_200_block1_GAM,             plot_PAA_200_block1_GAM,
             plot_PA_200_block2_GAM,     plot_PBG_200_block2_GAM,             plot_PAA_200_block2_GAM,
             plot_PA_200_clusters_GAM,   plot_PBG_200_clusters_GAM,           plot_PAA_200_clusters_GAM,
             plot_PA_200_KNNDM_GAM,      plot_PBG_200_KNNDM_GAM,              plot_PAA_200_KNNDM_GAM,
             plot_PA_200_random_GAM,     plot_PBG_200_random_GAM,             plot_PAA_200_random_GAM,
             plot_PA_200_block1_Lasso,   plot_PBG_200_block1_Lasso,           plot_PAA_200_block1_Lasso,
             plot_PA_200_block2_Lasso,   plot_PBG_200_block2_Lasso,           plot_PAA_200_block2_Lasso,
             plot_PA_200_clusters_Lasso, plot_PBG_200_clusters_Lasso,         plot_PAA_200_clusters_Lasso,
             plot_PA_200_KNNDM_Lasso,    plot_PBG_200_KNNDM_Lasso,            plot_PAA_200_KNNDM_Lasso,
             plot_PA_200_random_Lasso,   plot_PBG_200_random_Lasso,           plot_PAA_200_random_Lasso,
             plot_PA_200_block1_Maxent,  plot_PBG_200_block1_Maxent,          plot_PAA_200_block1_Maxent,
             plot_PA_200_block2_Maxent,  plot_PBG_200_block2_Maxent,          plot_PAA_200_block2_Maxent,
             plot_PA_200_clusters_Maxent,plot_PBG_200_clusters_Maxent,        plot_PAA_200_clusters_Maxent,
             plot_PA_200_KNNDM_Maxent,   plot_PBG_200_KNNDM_Maxent,           plot_PAA_200_KNNDM_Maxent,
             plot_PA_200_random_Maxent,  plot_PBG_200_random_Maxent,          plot_PAA_200_random_Maxent,
             plot_PA_200_block1_RF,      plot_PBG_200_block1_RF,              plot_PAA_200_block1_RF,
             plot_PA_200_block2_RF,      plot_PBG_200_block2_RF,              plot_PAA_200_block2_RF,
             plot_PA_200_clusters_RF,    plot_PBG_200_clusters_RF,            plot_PAA_200_clusters_RF,
             plot_PA_200_KNNDM_RF,       plot_PBG_200_KNNDM_RF,               plot_PAA_200_KNNDM_RF,
             plot_PA_200_random_RF,      plot_PBG_200_random_RF,              plot_PAA_200_random_RF,
             plot_PA_400_block1_BRT,     plot_PBG_400_block1_BRT,             plot_PAA_400_block1_BRT,
             plot_PA_400_block2_BRT,     plot_PBG_400_block2_BRT,             plot_PAA_400_block2_BRT,
             plot_PA_400_clusters_BRT,   plot_PBG_400_clusters_BRT,           plot_PAA_400_clusters_BRT,
             plot_PA_400_KNNDM_BRT,      plot_PBG_400_KNNDM_BRT,              plot_PAA_400_KNNDM_BRT,
             plot_PA_400_random_BRT,     plot_PBG_400_random_BRT,             plot_PAA_400_random_BRT,
             plot_PA_400_block1_GAM,     plot_PBG_400_block1_GAM,             plot_PAA_400_block1_GAM,
             plot_PA_400_block2_GAM,     plot_PBG_400_block2_GAM,             plot_PAA_400_block2_GAM,
             plot_PA_400_clusters_GAM,   plot_PBG_400_clusters_GAM,           plot_PAA_400_clusters_GAM,
             plot_PA_400_KNNDM_GAM,      plot_PBG_400_KNNDM_GAM,              plot_PAA_400_KNNDM_GAM,
             plot_PA_400_random_GAM,     plot_PBG_400_random_GAM,             plot_PAA_400_random_GAM,
             plot_PA_400_block1_Lasso,   plot_PBG_400_block1_Lasso,           plot_PAA_400_block1_Lasso,
             plot_PA_400_block2_Lasso,   plot_PBG_400_block2_Lasso,           plot_PAA_400_block2_Lasso,
             plot_PA_400_clusters_Lasso, plot_PBG_400_clusters_Lasso,         plot_PAA_400_clusters_Lasso,
             plot_PA_400_KNNDM_Lasso,    plot_PBG_400_KNNDM_Lasso,            plot_PAA_400_KNNDM_Lasso,
             plot_PA_400_random_Lasso,   plot_PBG_400_random_Lasso,           plot_PAA_400_random_Lasso,
             plot_PA_400_block1_Maxent,  plot_PBG_400_block1_Maxent,          plot_PAA_400_block1_Maxent,
             plot_PA_400_block2_Maxent,  plot_PBG_400_block2_Maxent,          plot_PAA_400_block2_Maxent,
             plot_PA_400_clusters_Maxent,plot_PBG_400_clusters_Maxent,        plot_PAA_400_clusters_Maxent,
             plot_PA_400_KNNDM_Maxent,   plot_PBG_400_KNNDM_Maxent,           plot_PAA_400_KNNDM_Maxent,
             plot_PA_400_random_Maxent,  plot_PBG_400_random_Maxent,          plot_PAA_400_random_Maxent,
             plot_PA_400_block1_RF,      plot_PBG_400_block1_RF,              plot_PAA_400_block1_RF,
             plot_PA_400_block2_RF,      plot_PBG_400_block2_RF,              plot_PAA_400_block2_RF,
             plot_PA_400_clusters_RF,    plot_PBG_400_clusters_RF,            plot_PAA_400_clusters_RF,
             plot_PA_400_KNNDM_RF,       plot_PBG_400_KNNDM_RF,               plot_PAA_400_KNNDM_RF,
             plot_PA_400_random_RF,      plot_PBG_400_random_RF,              plot_PAA_400_random_RF)
  
  rm(plot_PA_40_block1_BRT,      plot_PBG_40_block1_BRT,              plot_PAA_40_block1_BRT,
     plot_PA_40_block2_BRT,      plot_PBG_40_block2_BRT,              plot_PAA_40_block2_BRT,
     plot_PA_40_clusters_BRT,    plot_PBG_40_clusters_BRT,            plot_PAA_40_clusters_BRT,
     plot_PA_40_KNNDM_BRT,       plot_PBG_40_KNNDM_BRT,               plot_PAA_40_KNNDM_BRT,
     plot_PA_40_random_BRT,      plot_PBG_40_random_BRT,              plot_PAA_40_random_BRT,
     plot_PA_40_block1_GAM,      plot_PBG_40_block1_GAM,              plot_PAA_40_block1_GAM,
     plot_PA_40_block2_GAM,      plot_PBG_40_block2_GAM,              plot_PAA_40_block2_GAM,
     plot_PA_40_clusters_GAM,    plot_PBG_40_clusters_GAM,            plot_PAA_40_clusters_GAM,
     plot_PA_40_KNNDM_GAM,       plot_PBG_40_KNNDM_GAM,               plot_PAA_40_KNNDM_GAM,
     plot_PA_40_random_GAM,      plot_PBG_40_random_GAM,              plot_PAA_40_random_GAM,
     plot_PA_40_block1_Lasso,    plot_PBG_40_block1_Lasso,            plot_PAA_40_block1_Lasso,
     plot_PA_40_block2_Lasso,    plot_PBG_40_block2_Lasso,            plot_PAA_40_block2_Lasso,
     plot_PA_40_clusters_Lasso,  plot_PBG_40_clusters_Lasso,          plot_PAA_40_clusters_Lasso,
     plot_PA_40_KNNDM_Lasso,     plot_PBG_40_KNNDM_Lasso,             plot_PAA_40_KNNDM_Lasso,
     plot_PA_40_random_Lasso,    plot_PBG_40_random_Lasso,            plot_PAA_40_random_Lasso,
     plot_PA_40_block1_Maxent,   plot_PBG_40_block1_Maxent,           plot_PAA_40_block1_Maxent,
     plot_PA_40_block2_Maxent,   plot_PBG_40_block2_Maxent,           plot_PAA_40_block2_Maxent,
     plot_PA_40_clusters_Maxent, plot_PBG_40_clusters_Maxent,         plot_PAA_40_clusters_Maxent,
     plot_PA_40_KNNDM_Maxent,    plot_PBG_40_KNNDM_Maxent,            plot_PAA_40_KNNDM_Maxent,
     plot_PA_40_random_Maxent,   plot_PBG_40_random_Maxent,           plot_PAA_40_random_Maxent,
     plot_PA_40_block1_RF,       plot_PBG_40_block1_RF,               plot_PAA_40_block1_RF,
     plot_PA_40_block2_RF,       plot_PBG_40_block2_RF,               plot_PAA_40_block2_RF,
     plot_PA_40_clusters_RF,     plot_PBG_40_clusters_RF,             plot_PAA_40_clusters_RF,
     plot_PA_40_KNNDM_RF,        plot_PBG_40_KNNDM_RF,                plot_PAA_40_KNNDM_RF,
     plot_PA_40_random_RF,       plot_PBG_40_random_RF,               plot_PAA_40_random_RF,
     plot_PA_80_block1_BRT,      plot_PBG_80_block1_BRT,              plot_PAA_80_block1_BRT,
     plot_PA_80_block2_BRT,      plot_PBG_80_block2_BRT,              plot_PAA_80_block2_BRT,
     plot_PA_80_clusters_BRT,    plot_PBG_80_clusters_BRT,            plot_PAA_80_clusters_BRT,
     plot_PA_80_KNNDM_BRT,       plot_PBG_80_KNNDM_BRT,               plot_PAA_80_KNNDM_BRT,
     plot_PA_80_random_BRT,      plot_PBG_80_random_BRT,              plot_PAA_80_random_BRT,
     plot_PA_80_block1_GAM,      plot_PBG_80_block1_GAM,              plot_PAA_80_block1_GAM,
     plot_PA_80_block2_GAM,      plot_PBG_80_block2_GAM,              plot_PAA_80_block2_GAM,
     plot_PA_80_clusters_GAM,    plot_PBG_80_clusters_GAM,            plot_PAA_80_clusters_GAM,
     plot_PA_80_KNNDM_GAM,       plot_PBG_80_KNNDM_GAM,               plot_PAA_80_KNNDM_GAM,
     plot_PA_80_random_GAM,      plot_PBG_80_random_GAM,              plot_PAA_80_random_GAM,
     plot_PA_80_block1_Lasso,    plot_PBG_80_block1_Lasso,            plot_PAA_80_block1_Lasso,
     plot_PA_80_block2_Lasso,    plot_PBG_80_block2_Lasso,            plot_PAA_80_block2_Lasso,
     plot_PA_80_clusters_Lasso,  plot_PBG_80_clusters_Lasso,          plot_PAA_80_clusters_Lasso,
     plot_PA_80_KNNDM_Lasso,     plot_PBG_80_KNNDM_Lasso,             plot_PAA_80_KNNDM_Lasso,
     plot_PA_80_random_Lasso,    plot_PBG_80_random_Lasso,            plot_PAA_80_random_Lasso,
     plot_PA_80_block1_Maxent,   plot_PBG_80_block1_Maxent,           plot_PAA_80_block1_Maxent,
     plot_PA_80_block2_Maxent,   plot_PBG_80_block2_Maxent,           plot_PAA_80_block2_Maxent,
     plot_PA_80_clusters_Maxent, plot_PBG_80_clusters_Maxent,         plot_PAA_80_clusters_Maxent,
     plot_PA_80_KNNDM_Maxent,    plot_PBG_80_KNNDM_Maxent,            plot_PAA_80_KNNDM_Maxent,
     plot_PA_80_random_Maxent,   plot_PBG_80_random_Maxent,           plot_PAA_80_random_Maxent,
     plot_PA_80_block1_RF,       plot_PBG_80_block1_RF,               plot_PAA_80_block1_RF,
     plot_PA_80_block2_RF,       plot_PBG_80_block2_RF,               plot_PAA_80_block2_RF,
     plot_PA_80_clusters_RF,     plot_PBG_80_clusters_RF,             plot_PAA_80_clusters_RF,
     plot_PA_80_KNNDM_RF,        plot_PBG_80_KNNDM_RF,                plot_PAA_80_KNNDM_RF,
     plot_PA_80_random_RF,       plot_PBG_80_random_RF,               plot_PAA_80_random_RF,
     plot_PA_120_block1_BRT,     plot_PBG_120_block1_BRT,             plot_PAA_120_block1_BRT,
     plot_PA_120_block2_BRT,     plot_PBG_120_block2_BRT,             plot_PAA_120_block2_BRT,
     plot_PA_120_clusters_BRT,   plot_PBG_120_clusters_BRT,           plot_PAA_120_clusters_BRT,
     plot_PA_120_KNNDM_BRT,      plot_PBG_120_KNNDM_BRT,              plot_PAA_120_KNNDM_BRT,
     plot_PA_120_random_BRT,     plot_PBG_120_random_BRT,             plot_PAA_120_random_BRT,
     plot_PA_120_block1_GAM,     plot_PBG_120_block1_GAM,             plot_PAA_120_block1_GAM,
     plot_PA_120_block2_GAM,     plot_PBG_120_block2_GAM,             plot_PAA_120_block2_GAM,
     plot_PA_120_clusters_GAM,   plot_PBG_120_clusters_GAM,           plot_PAA_120_clusters_GAM,
     plot_PA_120_KNNDM_GAM,      plot_PBG_120_KNNDM_GAM,              plot_PAA_120_KNNDM_GAM,
     plot_PA_120_random_GAM,     plot_PBG_120_random_GAM,             plot_PAA_120_random_GAM,
     plot_PA_120_block1_Lasso,   plot_PBG_120_block1_Lasso,           plot_PAA_120_block1_Lasso,
     plot_PA_120_block2_Lasso,   plot_PBG_120_block2_Lasso,           plot_PAA_120_block2_Lasso,
     plot_PA_120_clusters_Lasso, plot_PBG_120_clusters_Lasso,         plot_PAA_120_clusters_Lasso,
     plot_PA_120_KNNDM_Lasso,    plot_PBG_120_KNNDM_Lasso,            plot_PAA_120_KNNDM_Lasso,
     plot_PA_120_random_Lasso,   plot_PBG_120_random_Lasso,           plot_PAA_120_random_Lasso,
     plot_PA_120_block1_Maxent,  plot_PBG_120_block1_Maxent,          plot_PAA_120_block1_Maxent,
     plot_PA_120_block2_Maxent,  plot_PBG_120_block2_Maxent,          plot_PAA_120_block2_Maxent,
     plot_PA_120_clusters_Maxent,plot_PBG_120_clusters_Maxent,        plot_PAA_120_clusters_Maxent,
     plot_PA_120_KNNDM_Maxent,   plot_PBG_120_KNNDM_Maxent,           plot_PAA_120_KNNDM_Maxent,
     plot_PA_120_random_Maxent,  plot_PBG_120_random_Maxent,          plot_PAA_120_random_Maxent,
     plot_PA_120_block1_RF,      plot_PBG_120_block1_RF,              plot_PAA_120_block1_RF,
     plot_PA_120_block2_RF,      plot_PBG_120_block2_RF,              plot_PAA_120_block2_RF,
     plot_PA_120_clusters_RF,    plot_PBG_120_clusters_RF,            plot_PAA_120_clusters_RF,
     plot_PA_120_KNNDM_RF,       plot_PBG_120_KNNDM_RF,               plot_PAA_120_KNNDM_RF,
     plot_PA_120_random_RF,      plot_PBG_120_random_RF,              plot_PAA_120_random_RF,
     plot_PA_160_block1_BRT,     plot_PBG_160_block1_BRT,             plot_PAA_160_block1_BRT,
     plot_PA_160_block2_BRT,     plot_PBG_160_block2_BRT,             plot_PAA_160_block2_BRT,
     plot_PA_160_clusters_BRT,   plot_PBG_160_clusters_BRT,           plot_PAA_160_clusters_BRT,
     plot_PA_160_KNNDM_BRT,      plot_PBG_160_KNNDM_BRT,              plot_PAA_160_KNNDM_BRT,
     plot_PA_160_random_BRT,     plot_PBG_160_random_BRT,             plot_PAA_160_random_BRT,
     plot_PA_160_block1_GAM,     plot_PBG_160_block1_GAM,             plot_PAA_160_block1_GAM,
     plot_PA_160_block2_GAM,     plot_PBG_160_block2_GAM,             plot_PAA_160_block2_GAM,
     plot_PA_160_clusters_GAM,   plot_PBG_160_clusters_GAM,           plot_PAA_160_clusters_GAM,
     plot_PA_160_KNNDM_GAM,      plot_PBG_160_KNNDM_GAM,              plot_PAA_160_KNNDM_GAM,
     plot_PA_160_random_GAM,     plot_PBG_160_random_GAM,             plot_PAA_160_random_GAM,
     plot_PA_160_block1_Lasso,   plot_PBG_160_block1_Lasso,           plot_PAA_160_block1_Lasso,
     plot_PA_160_block2_Lasso,   plot_PBG_160_block2_Lasso,           plot_PAA_160_block2_Lasso,
     plot_PA_160_clusters_Lasso, plot_PBG_160_clusters_Lasso,         plot_PAA_160_clusters_Lasso,
     plot_PA_160_KNNDM_Lasso,    plot_PBG_160_KNNDM_Lasso,            plot_PAA_160_KNNDM_Lasso,
     plot_PA_160_random_Lasso,   plot_PBG_160_random_Lasso,           plot_PAA_160_random_Lasso,
     plot_PA_160_block1_Maxent,  plot_PBG_160_block1_Maxent,          plot_PAA_160_block1_Maxent,
     plot_PA_160_block2_Maxent,  plot_PBG_160_block2_Maxent,          plot_PAA_160_block2_Maxent,
     plot_PA_160_clusters_Maxent,plot_PBG_160_clusters_Maxent,        plot_PAA_160_clusters_Maxent,
     plot_PA_160_KNNDM_Maxent,   plot_PBG_160_KNNDM_Maxent,           plot_PAA_160_KNNDM_Maxent,
     plot_PA_160_random_Maxent,  plot_PBG_160_random_Maxent,          plot_PAA_160_random_Maxent,
     plot_PA_160_block1_RF,      plot_PBG_160_block1_RF,              plot_PAA_160_block1_RF,
     plot_PA_160_block2_RF,      plot_PBG_160_block2_RF,              plot_PAA_160_block2_RF,
     plot_PA_160_clusters_RF,    plot_PBG_160_clusters_RF,            plot_PAA_160_clusters_RF,
     plot_PA_160_KNNDM_RF,       plot_PBG_160_KNNDM_RF,               plot_PAA_160_KNNDM_RF,
     plot_PA_160_random_RF,      plot_PBG_160_random_RF,              plot_PAA_160_random_RF,
     plot_PA_200_block1_BRT,     plot_PBG_200_block1_BRT,             plot_PAA_200_block1_BRT,
     plot_PA_200_block2_BRT,     plot_PBG_200_block2_BRT,             plot_PAA_200_block2_BRT,
     plot_PA_200_clusters_BRT,   plot_PBG_200_clusters_BRT,           plot_PAA_200_clusters_BRT,
     plot_PA_200_KNNDM_BRT,      plot_PBG_200_KNNDM_BRT,              plot_PAA_200_KNNDM_BRT,
     plot_PA_200_random_BRT,     plot_PBG_200_random_BRT,             plot_PAA_200_random_BRT,
     plot_PA_200_block1_GAM,     plot_PBG_200_block1_GAM,             plot_PAA_200_block1_GAM,
     plot_PA_200_block2_GAM,     plot_PBG_200_block2_GAM,             plot_PAA_200_block2_GAM,
     plot_PA_200_clusters_GAM,   plot_PBG_200_clusters_GAM,           plot_PAA_200_clusters_GAM,
     plot_PA_200_KNNDM_GAM,      plot_PBG_200_KNNDM_GAM,              plot_PAA_200_KNNDM_GAM,
     plot_PA_200_random_GAM,     plot_PBG_200_random_GAM,             plot_PAA_200_random_GAM,
     plot_PA_200_block1_Lasso,   plot_PBG_200_block1_Lasso,           plot_PAA_200_block1_Lasso,
     plot_PA_200_block2_Lasso,   plot_PBG_200_block2_Lasso,           plot_PAA_200_block2_Lasso,
     plot_PA_200_clusters_Lasso, plot_PBG_200_clusters_Lasso,         plot_PAA_200_clusters_Lasso,
     plot_PA_200_KNNDM_Lasso,    plot_PBG_200_KNNDM_Lasso,            plot_PAA_200_KNNDM_Lasso,
     plot_PA_200_random_Lasso,   plot_PBG_200_random_Lasso,           plot_PAA_200_random_Lasso,
     plot_PA_200_block1_Maxent,  plot_PBG_200_block1_Maxent,          plot_PAA_200_block1_Maxent,
     plot_PA_200_block2_Maxent,  plot_PBG_200_block2_Maxent,          plot_PAA_200_block2_Maxent,
     plot_PA_200_clusters_Maxent,plot_PBG_200_clusters_Maxent,        plot_PAA_200_clusters_Maxent,
     plot_PA_200_KNNDM_Maxent,   plot_PBG_200_KNNDM_Maxent,           plot_PAA_200_KNNDM_Maxent,
     plot_PA_200_random_Maxent,  plot_PBG_200_random_Maxent,          plot_PAA_200_random_Maxent,
     plot_PA_200_block1_RF,      plot_PBG_200_block1_RF,              plot_PAA_200_block1_RF,
     plot_PA_200_block2_RF,      plot_PBG_200_block2_RF,              plot_PAA_200_block2_RF,
     plot_PA_200_clusters_RF,    plot_PBG_200_clusters_RF,            plot_PAA_200_clusters_RF,
     plot_PA_200_KNNDM_RF,       plot_PBG_200_KNNDM_RF,               plot_PAA_200_KNNDM_RF,
     plot_PA_200_random_RF,      plot_PBG_200_random_RF,              plot_PAA_200_random_RF,
     plot_PA_400_block1_BRT,     plot_PBG_400_block1_BRT,             plot_PAA_400_block1_BRT,
     plot_PA_400_block2_BRT,     plot_PBG_400_block2_BRT,             plot_PAA_400_block2_BRT,
     plot_PA_400_clusters_BRT,   plot_PBG_400_clusters_BRT,           plot_PAA_400_clusters_BRT,
     plot_PA_400_KNNDM_BRT,      plot_PBG_400_KNNDM_BRT,              plot_PAA_400_KNNDM_BRT,
     plot_PA_400_random_BRT,     plot_PBG_400_random_BRT,             plot_PAA_400_random_BRT,
     plot_PA_400_block1_GAM,     plot_PBG_400_block1_GAM,             plot_PAA_400_block1_GAM,
     plot_PA_400_block2_GAM,     plot_PBG_400_block2_GAM,             plot_PAA_400_block2_GAM,
     plot_PA_400_clusters_GAM,   plot_PBG_400_clusters_GAM,           plot_PAA_400_clusters_GAM,
     plot_PA_400_KNNDM_GAM,      plot_PBG_400_KNNDM_GAM,              plot_PAA_400_KNNDM_GAM,
     plot_PA_400_random_GAM,     plot_PBG_400_random_GAM,             plot_PAA_400_random_GAM,
     plot_PA_400_block1_Lasso,   plot_PBG_400_block1_Lasso,           plot_PAA_400_block1_Lasso,
     plot_PA_400_block2_Lasso,   plot_PBG_400_block2_Lasso,           plot_PAA_400_block2_Lasso,
     plot_PA_400_clusters_Lasso, plot_PBG_400_clusters_Lasso,         plot_PAA_400_clusters_Lasso,
     plot_PA_400_KNNDM_Lasso,    plot_PBG_400_KNNDM_Lasso,            plot_PAA_400_KNNDM_Lasso,
     plot_PA_400_random_Lasso,   plot_PBG_400_random_Lasso,           plot_PAA_400_random_Lasso,
     plot_PA_400_block1_Maxent,  plot_PBG_400_block1_Maxent,          plot_PAA_400_block1_Maxent,
     plot_PA_400_block2_Maxent,  plot_PBG_400_block2_Maxent,          plot_PAA_400_block2_Maxent,
     plot_PA_400_clusters_Maxent,plot_PBG_400_clusters_Maxent,        plot_PAA_400_clusters_Maxent,
     plot_PA_400_KNNDM_Maxent,   plot_PBG_400_KNNDM_Maxent,           plot_PAA_400_KNNDM_Maxent,
     plot_PA_400_random_Maxent,  plot_PBG_400_random_Maxent,          plot_PAA_400_random_Maxent,
     plot_PA_400_block1_RF,      plot_PBG_400_block1_RF,              plot_PAA_400_block1_RF,
     plot_PA_400_block2_RF,      plot_PBG_400_block2_RF,              plot_PAA_400_block2_RF,
     plot_PA_400_clusters_RF,    plot_PBG_400_clusters_RF,            plot_PAA_400_clusters_RF,
     plot_PA_400_KNNDM_RF,       plot_PBG_400_KNNDM_RF,               plot_PAA_400_KNNDM_RF,
     plot_PA_400_random_RF,      plot_PBG_400_random_RF,              plot_PAA_400_random_RF)
  
  p1=egg::ggarrange(plots=plots[1:15],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p1, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",40,".png"), dpi = 300, width = 16, height = 25)
  
  p2=egg::ggarrange(plots=plots[16:30],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p2, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",40,".png"), dpi = 300, width = 16, height = 25)
  
  p3=egg::ggarrange(plots=plots[31:45],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p3, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",40,".png"), dpi = 300, width = 16, height = 25)
  
  p4=egg::ggarrange(plots=plots[46:60],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p4, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",40,".png"), dpi = 300, width = 16, height = 25)
  
  
  p5=egg::ggarrange(plots=plots[61:75],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p5, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",40,".png"), dpi = 300, width = 16, height = 25)
  
  p6=egg::ggarrange(plots=plots[76:90],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p6, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",80,".png"), dpi = 300, width = 16, height = 25)
  
  p7=egg::ggarrange(plots=plots[91:105],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p7, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",80,".png"), dpi = 300, width = 16, height = 25)
  
  p8=egg::ggarrange(plots=plots[106:120],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p8, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",80,".png"), dpi = 300, width = 16, height = 25)
  
  p9=egg::ggarrange(plots=plots[121:135],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p9, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",80,".png"), dpi = 300, width = 16, height = 25)
  
  
  p10=egg::ggarrange(plots=plots[136:150],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p10, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",80,".png"), dpi = 300, width = 16, height = 25)
  
  p11=egg::ggarrange(plots=plots[151:165],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p11, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",120,".png"), dpi = 300, width = 16, height = 25)
  
  p12=egg::ggarrange(plots=plots[166:180],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p12, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",120,".png"), dpi = 300, width = 16, height = 25)
  
  p13=egg::ggarrange(plots=plots[181:195],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p13, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",120,".png"), dpi = 300, width = 16, height = 25)
  
  p14=egg::ggarrange(plots=plots[196:210],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p14, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",120,".png"), dpi = 300, width = 16, height = 25)
  
  
  p15=egg::ggarrange(plots=plots[211:225],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p15, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",120,".png"), dpi = 300, width = 16, height = 25)
  
  p16=egg::ggarrange(plots=plots[226:240],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p16, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",160,".png"), dpi = 300, width = 16, height = 25)
  
  p17=egg::ggarrange(plots=plots[241:255],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p17, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",160,".png"), dpi = 300, width = 16, height = 25)
  
  p18=egg::ggarrange(plots=plots[256:270],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p18, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",160,".png"), dpi = 300, width = 16, height = 25)
  
  p19=egg::ggarrange(plots=plots[271:285],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p19, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",160,".png"), dpi = 300, width = 16, height = 25)
  
  
  p20=egg::ggarrange(plots=plots[286:300],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p20, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",160,".png"), dpi = 300, width = 16, height = 25)
  
  p21=egg::ggarrange(plots=plots[301:315],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p21, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",200,".png"), dpi = 300, width = 16, height = 25)
  
  p22=egg::ggarrange(plots=plots[316:330],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p22, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",200,".png"), dpi = 300, width = 16, height = 25)
  
  p23=egg::ggarrange(plots=plots[331:345],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p23, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",200,".png"), dpi = 300, width = 16, height = 25)
  
  p24=egg::ggarrange(plots=plots[346:360],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p24, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",200,".png"), dpi = 300, width = 16, height = 25)
  
  
  p25=egg::ggarrange(plots=plots[361:375],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p25, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",200,".png"), dpi = 300, width = 16, height = 25)
  
  p26=egg::ggarrange(plots=plots[376:390],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p26, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_BRT_size",400,".png"), dpi = 300, width = 16, height = 25)
  
  p27=egg::ggarrange(plots=plots[391:405],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p27, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_GAM_size",400,".png"), dpi = 300, width = 16, height = 25)
  
  p28=egg::ggarrange(plots=plots[406:420],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p28, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Lasso_size",400,".png"), dpi = 300, width = 16, height = 25)
  
  p29=egg::ggarrange(plots=plots[421:435],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p29, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_Maxent_size",400,".png"), dpi = 300, width = 16, height = 25)
  
  
  p30=egg::ggarrange(plots=plots[436:450],ncol=3,nrow=5, left=textGrob(
    "Pearson correlation between probability of occurrence and artificial distribution map",
    gp = gpar(fontsize = 24),
    rot = 90))
  ggsave(p30, filename = paste0("images/",nameRun,"/resultPlotsSizeModelPoints/", m, "_RF_size",400,".png"), dpi = 300, width = 16, height = 25)
  
  
  
  rm(plots,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30)
}
