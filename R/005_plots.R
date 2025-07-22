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

# 2 - load data ####
#------------------#

data=readRDS("data/run2/results.RDS")

data <- data %>% filter(is.finite(metric))%>%
  mutate(trueCor = ifelse(trueCor < 0, 0, trueCor))

#index calculation

data=data %>%
  dplyr::mutate(
    metric = rowMeans(cbind(Spec, COR, PCC, 1 - MAE), na.rm = TRUE)
  )

if(!dir.exists("images/resultPlots")) dir.create("images/resultPlots")

# 3 - plot by method ####
#-----------------------#

metric_names <- c("metric", "AUC", "COR", "Spec", "Sens", "Kappa", 
                  "PCC", "TSS", "PRG", "MAE", "BIAS")

plot_list <- list()

for (m in metric_names) {
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", m))
  
  metrics_df <- data %>%
    group_by(method) %>%
    summarise(
      r_squared = summary(lm(formula, data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(formula, data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(formula, data = pick(everything())))[[m]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p <- ggplot(data, aes(x = .data[[m]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(method)) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  ggsave(p, filename = paste0("images/resultPlots/", m, "_byMethod.png"), dpi = 300, width = 16, height = 8)
  plot_list[[m]] <- p
  
  rm(p, metrics_df)
}

for (i in c(1,4,7,10)){
  
  filtered_plots <- Filter(Negate(is.null), plot_list[i:(i + 2)])
  p=gridExtra::grid.arrange(grobs = filtered_plots, ncol=1, nrow=3)
  ggsave(p, filename = paste0("images/resultPlots/byMethod",i,".png"), dpi = 300, width = 16, height = 16)
  rm(p)
};rm(plot_list,formula, filtered_plots,i)

# 4 - plot by size ####
#---------------------#

plot_list <- list()

for (m in metric_names) {
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", m))
  
  metrics_df <- data %>%
    group_by(method,size) %>%
    summarise(
      r_squared = summary(lm(formula, data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(formula, data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(formula, data = pick(everything())))[[m]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p=ggplot(data, aes(x = .data[[m]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(method, size)) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  ggsave(p, filename = paste0("images/resultPlots/", m, "_byMethodSize.png"), dpi = 300, width = 16, height = 16)
  plot_list[[m]] <- p
  
  rm(p, metrics_df)
};rm(plot_list,formula)


# 5 - plot by model ####
#----------------------#

plot_list <- list()

for (m in metric_names) {
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", m))
  
  metrics_df <- data %>%
    group_by(method,model) %>%
    summarise(
      r_squared = summary(lm(formula, data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(formula, data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(formula, data = pick(everything())))[[m]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p=ggplot(data, aes(x = .data[[m]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(method, model)) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  ggsave(p, filename = paste0("images/resultPlots/", m, "_byMethodModel.png"), dpi = 300, width = 16, height = 16)
  plot_list[[m]] <- p
  
  rm(p, metrics_df)
};rm(plot_list,formula)

# 6 - plot by points ####
#-----------------------#

plot_list <- list()

for (m in metric_names) {
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", m))
  
  metrics_df <- data %>%
    group_by(method,points) %>%
    summarise(
      r_squared = summary(lm(formula, data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(formula, data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(formula, data = pick(everything())))[[m]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p=ggplot(data, aes(x = .data[[m]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(method, points)) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  ggsave(p, filename = paste0("images/resultPlots/", m, "_byMethodPoints.png"), dpi = 300, width = 16, height = 16)
  plot_list[[m]] <- p
  
  rm(p, metrics_df)
};rm(plot_list,formula)



############################################



data$testData_bin <- cut(
  data$noPresencePoints,
  breaks = seq(0, max(data$noPresencePoints, na.rm = TRUE) + 10, by = 20),
  right = TRUE,    # include the upper bound in the interval
  include.lowest = TRUE,
  labels = paste(seq(1, max(data$noPresencePoints, na.rm = TRUE), by = 20),
                 seq(10, max(data$noPresencePoints, na.rm = TRUE) + 10, by = 20),
                 sep = "–")
)


