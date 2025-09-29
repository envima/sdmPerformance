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
library(lmodel2)


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

if(!dir.exists(paste0("images/",nameRun,"/resultPlots"))) dir.create(paste0("images/",nameRun,"/resultPlots"), recursive=T)

# 3 - plot with reduced major axis regression ####
#------------------------------------------------#

#reduced major axis same as standard major axis (SMA)
# http://stratigrafia.org/8370/lecturenotes/regression.html
# https://stackoverflow.com/questions/49196327/is-there-a-difference-between-reduced-major-axis-regression-and-ranged-major-axi

metric_names <- c( "AUC_scaled", "COR_scaled", "Spec_scaled", "Sens_scaled", "Kappa_scaled", 
                   "PCC_scaled", "TSS_scaled", "PRG_scaled", "SEDI_scaled", "SBI_m_scaled")

for (m in metric_names) {
  
  for(method_plot in c("PAA", "PA", "PBG")){
    
    data2=data%>%dplyr::filter(method == method_plot)
    # Fit RMA regression for this metric
    rma_fit <- lmodel2::lmodel2(
      formula = as.formula(paste("trueCor_scaled ~", m)),
      data = data2,
      nperm = 0 # no permutation test
    )
    
    # calculate slope and intercept     
    slope_rma <- rma_fit$regression.results$Slope[rma_fit$regression.results$Method=="SMA"]
    intercept_rma <- rma_fit$regression.results$Intercept[rma_fit$regression.results$Method=="SMA"]
    
    
    
    # Compute predicted values along the RMA line
    data2$y_rma <- intercept_rma + slope_rma * data2[[m]]
    
    # Metrics relative to RMA
    metrics_df <- data2 %>%
      # group_by(method) %>%
      summarise(
        rmse_rma = sqrt(mean((trueCor_scaled - y_rma)^2, na.rm = TRUE)),
        mae_rma  = mean(abs(trueCor_scaled - y_rma), na.rm = TRUE),
        r_squared_rma = 1 - sum((trueCor_scaled - y_rma)^2, na.rm = TRUE) / sum((trueCor_scaled - mean(trueCor_scaled, na.rm = TRUE))^2, na.rm = TRUE),
        # mean_intercept_dev_rma = mean(trueCor_scaled - y_rma, na.rm = TRUE),
        xpos = 0,
        ypos = 0.95,
        .groups = "drop"
      )
    
    metrics_df$slope_dev <- slope_rma-1
    metrics_df$intercept_dev <- intercept_rma-0
    
    p <- ggplot(data2, aes(x = .data[[m]], y = trueCor_scaled)) +
      geom_point(size = 0.3, colour="cornflowerblue") +
      geom_abline(slope = 1, intercept = 0, color = "deeppink3", linetype = "dashed", linewidth=1) + # Winkelhalbierende
      geom_abline(slope = slope_rma, intercept = intercept_rma, linewidth=1) + # RMA line
      facet_wrap(vars(method)) +
      geom_text(
        data = metrics_df,
        aes(x = xpos, y = ypos,
            label = paste0("R² = ", round(r_squared_rma, 2),
                           "\nRMSE = ", round(rmse_rma, 2),
                           "\nMAE = ", round(mae_rma, 2),
                           "\nΔslope = ", round(slope_dev, 2),
                           "\nΔintercept = ", round(intercept_dev, 2)
            )),
        inherit.aes = FALSE,
        size = 5,
        hjust = 0
      ) +
      xlim(0, 1) + ylim(-0.01, 1) +
      xlab(m) +
      ylab("Pearson correlation between suitability \nraster and artificial distribution map") +
      theme_minimal(base_size = 15)+theme(legend.position="none")
    
    assign(paste0("plot_",method_plot), p)
    rm(data2, rma_fit,intercept_rma, slope_rma,p,metrics_df, slope_dev,intercept_dev)
    
  }
  
  
  plots=list(plot_PA, plot_PBG, plot_PAA)
  remove_axis <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  plots[-c(1,4)] <- lapply(plots[-c(1,4)] , function(.p) .p + remove_axis)
  
  p=egg::ggarrange(plots=plots,ncol=3)
  ggsave(p, filename = paste0("images/",nameRun,"/resultPlots/", m, "_byMethod.png"), dpi = 300, width = 16, height = 8)
  rm(plot_PA, plot_PBG, plot_PAA,plots,p,remove_axis)
}



