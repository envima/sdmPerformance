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
  setwd("M:/user/bald/SDM/sdmPerformance")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}



# create a unique name for different runs
nameRun <- paste0("run", 3)


source(paste0("R/",nameRun,"/functions/scale_metric.R"))

# 2 - load data ####
#------------------#


data=readRDS(paste0("data/",nameRun,"/results.RDS"))%>%dplyr::filter(noPresencePoints > 20)
data$points <-factor(data$points, levels = c(40, 80, 120, 160, 200, 400))

# List of metrics to scale
metric_names <- c("AUC", "COR", "Spec", "Sens", "Kappa", 
                  "PCC", "TSS", "PRG", "MAE", "trueCor")

# Apply scaling function to each metric and create new "_scaled" columns
data <- data %>%
  mutate(across(all_of(metric_names),
                .fns = ~ scale_metric(., cur_column()),
                .names = "{.col}_scaled"))


data=data %>%
  dplyr::mutate(
    metric_scaled = rowMeans(cbind( COR_scaled, TSS_scaled, MAE_scaled), na.rm = TRUE)
  )


if(!dir.exists(paste0("images/",nameRun,"/resultPlots"))) dir.create(paste0("images/",nameRun,"/resultPlots"), recursive=T)

# 3 - plot with bisection by method ####
#--------------------------------------#

metric_names <- c( "metric_scaled","AUC_scaled", "COR_scaled", "Spec_scaled", "Sens_scaled", "Kappa_scaled", 
                   "PCC_scaled", "TSS_scaled", "PRG_scaled", "MAE_scaled")



for(j in c("points", "model", "size","species","replicate")){
  plot_list <- list()
  for (m in metric_names) {
    
    # Fit models dynamically using as.formula
    formula <- as.formula(paste("trueCor_scaled ~", m))
    # Berechnung der Kennzahlen bezogen auf die Winkelhalbierende y = x
    metrics_df <- data %>%
      group_by(method) %>%
      summarise(
        rmse = sqrt(mean((trueCor_scaled - .data[[m]])^2, na.rm = TRUE)),
        mae  = mean(abs(trueCor_scaled - .data[[m]]), na.rm = TRUE),
        r_squared = 1 - sum((trueCor_scaled - .data[[m]])^2, na.rm = TRUE) / sum((trueCor_scaled - mean(trueCor_scaled, na.rm = TRUE))^2, na.rm = TRUE),
        mean_intercept_deviation = mean(trueCor_scaled - .data[[m]], na.rm = TRUE),
        xpos = 0,
        ypos = 0.1,
        .groups = "drop"
      )
    
    p<-ggplot(data, aes(x = .data[[m]], y = trueCor_scaled, color=.data[[j]])) +
      geom_point(size=0.01)+
      # geom_hex(bins = 140, alpha = 1) +
      #scale_fill_viridis_c(name = "Point density") +
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") + # Winkelhalbierende
      facet_wrap(vars(method)) +
      geom_text(
        data = metrics_df,
        aes(x = xpos, y = ypos,
            label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2),  "\nMAE = ", round(mae, 2), 
                           "\nMean Intercept Dev. = ", round(mean_intercept_deviation, 2))),
        inherit.aes = FALSE,
        size = 5,
        hjust = 0
      ) +
      xlim(0, 1) + ylim(-0.01, 1) +
      xlab(m) +
      ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
      theme_minimal(base_size = 15)
    
    
    
    
    ggsave(p, filename = paste0("images/",nameRun,"/resultPlots/", m, "_byMethod_color",j,".png"), dpi = 300, width = 16, height = 8)
    plot_list[[m]] <- p
    
    rm(p, metrics_df)
  }
}


# 4 - plot with reduced major axis regression ####
#------------------------------------------------#

#reduced major axis same as standard major axis (SMA)
# http://stratigrafia.org/8370/lecturenotes/regression.html
# https://stackoverflow.com/questions/49196327/is-there-a-difference-between-reduced-major-axis-regression-and-ranged-major-axi


metric_names <- c( "metric_scaled","AUC_scaled", "COR_scaled", "Spec_scaled", "Sens_scaled", "Kappa_scaled", 
                   "PCC_scaled", "TSS_scaled", "PRG_scaled", "MAE_scaled")


for(j in c("points", "model", "size","species","replicate")){
  plot_list <- list()
  for (m in metric_names) {
    
    
    # Fit RMA regression for this metric
    rma_fit_paa <- lmodel2::lmodel2(
      formula = as.formula(paste("trueCor_scaled ~", m)),
      data = data,
      
      nperm = 0 # no permutation test
    )
    
    metrics_df <- data %>%
      group_by(method) %>%
      summarise(
        rma_fit = list(lmodel2::lmodel2(as.formula(paste("trueCor_scaled ~", m)),
                                        data = cur_data(), nperm = 0)),
        .groups = "drop"
      )
    
    slope_rma <- rma_fit$regression.results$Slope[rma_fit$regression.results$Method=="SMA"]
    intercept_rma <- rma_fit$regression.results$Intercept[rma_fit$regression.results$Method=="SMA"]
    
    # Compute predicted values along the RMA line
    data$y_rma <- intercept_rma + slope_rma * data[[m]]
    
    # Metrics relative to RMA
    metrics_df <- data %>%
      group_by(method) %>%
      summarise(
        rmse_rma = sqrt(mean((trueCor_scaled - y_rma)^2, na.rm = TRUE)),
        mae_rma  = mean(abs(trueCor_scaled - y_rma), na.rm = TRUE),
        r_squared_rma = 1 - sum((trueCor_scaled - y_rma)^2, na.rm = TRUE) / sum((trueCor_scaled - mean(trueCor_scaled, na.rm = TRUE))^2, na.rm = TRUE),
        mean_intercept_dev_rma = mean(trueCor_scaled - y_rma, na.rm = TRUE),
        xpos = 0,
        ypos = 0.1,
        .groups = "drop"
      )
    
    
    p <- ggplot(data, aes(x = .data[[m]], y = trueCor_scaled, color=.data[[j]])) +
      geom_point(size = 0.01) +
      geom_abline(slope = slope_rma, intercept = intercept_rma, color = "blue", linetype = "dashed") + # RMA line
      facet_wrap(vars(method)) +
      geom_text(
        data = metrics_df,
        aes(x = xpos, y = ypos,
            label = paste0("R² = ", round(r_squared_rma, 2),
                           "\nRMSE = ", round(rmse_rma, 2),
                           "\nMAE = ", round(mae_rma, 2),
                           "\nMean Intercept Dev. = ", round(mean_intercept_dev_rma, 2))),
        inherit.aes = FALSE,
        size = 5,
        hjust = 0
      ) +
      xlim(0, 1) + ylim(-0.01, 1) +
      xlab(m) +
      ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
      theme_minimal(base_size = 15)
    
    
    
    
    ggsave(p, filename = paste0("images/",nameRun,"/resultPlots/", m, "_byMethod_color",j,".png"), dpi = 300, width = 16, height = 8)
    plot_list[[m]] <- p
    
    rm(p, metrics_df)
  }
}


