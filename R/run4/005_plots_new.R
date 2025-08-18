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

data=readRDS("data/run4/results.RDS")#%>%dplyr::filter(noPresencePoints > 20)

#data <- data %>% filter(is.finite(metric))%>%
#  mutate(trueCor = ifelse(trueCor < 0, 0, trueCor))


# Definiere Hilfsfunktionen


# List of metrics to scale
metric_names <- c("AUC", "COR", "Spec", "Sens", "Kappa", 
                  "PCC", "TSS", "PRG", "MAE", "BIAS")

# Apply scaling function to each metric and create new "_scaled" columns
data <- data %>%
  mutate(across(all_of(metric_names),
                .fns = ~ scale_metric(., cur_column()),
                .names = "{.col}_scaled"))

data=data%>%dplyr::mutate(trueCor_scaled=scale_metric(trueCor,"COR"))
#index calculation

data=data %>%
  dplyr::mutate(
    metric_scaled = rowMeans(cbind( COR_scaled, PCC_scaled, MAE_scaled), na.rm = TRUE)
  )

data$points <-factor(data$points, levels = c(40, 80, 120, 160, 200, 400))

if(!dir.exists("images/run4/resultPlots")) dir.create("images/run4/resultPlots", recursive=T)

# 3 - plot by method ####
#-----------------------#

metric_names <- c( "metric_scaled","AUC_scaled", "COR_scaled", "Spec_scaled", "Sens_scaled", "Kappa_scaled", 
                   "PCC_scaled", "TSS_scaled", "PRG_scaled", "MAE_scaled", "BIAS_scaled")

plot_list <- list()

for (m in metric_names) {
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", m))
  # Berechnung der Kennzahlen bezogen auf die Winkelhalbierende y = x
  metrics_df <- data %>%
    group_by(method) %>%
    summarise(
      rmse = sqrt(mean((trueCor - .data[[m]])^2, na.rm = TRUE)),
      r_squared = 1 - sum((trueCor - .data[[m]])^2, na.rm = TRUE) / sum((trueCor - mean(trueCor, na.rm = TRUE))^2, na.rm = TRUE),
      mean_intercept_deviation = mean(trueCor - .data[[m]], na.rm = TRUE),
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p<-ggplot(data, aes(x = .data[[m]], y = trueCor, color=points)) +
    geom_point(size=1)+
    # geom_hex(bins = 140, alpha = 1) +
    #scale_fill_viridis_c(name = "Point density") +
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") + # Winkelhalbierende
    facet_wrap(vars(method)) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos,
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), 
                         "\nMean Intercept Dev. = ", round(mean_intercept_deviation, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  
  
  
  ggsave(p, filename = paste0("images/run4/resultPlots/", m, "_byMethod_colorPoints.png"), dpi = 300, width = 16, height = 8)
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

library(paletteer)
library(forcats)

# 7 - plots on overall model performance ####
#-------------------------------------------#


data <- data %>%
  mutate(
    model = fct_reorder(model, trueCor, .fun = median, na.rm = TRUE),
    size = fct_reorder(size, trueCor, .fun = median, na.rm = TRUE),
    points = fct_reorder(points, trueCor, .fun = median, na.rm = TRUE)
  )

model<- ggplot(data, aes(x=model,y=trueCor, fill=model))+
  geom_boxplot()+ facet_wrap(vars(method))+scale_fill_paletteer_d("nationalparkcolors::Acadia", name = "Model \nalgorithm") +
  xlab("Model algorithm")+ylab("Pearson correlation between suitability \nraster and prediction map of virtual species")+
  theme_minimal(base_size = 15)


size<- ggplot(data, aes(x=size,y=trueCor, fill=size))+
  geom_boxplot()+ facet_wrap(vars(method))+scale_fill_paletteer_d("nationalparkcolors::Acadia", name = "Data \nseparation \nstrategy") +
  xlab("Data separation strategy")+ylab("Pearson correlation between suitability \nraster and prediction map of virtual species")+
  theme_minimal(base_size = 15)

points<- ggplot(data, aes(x=points,y=trueCor, fill=points))+
  geom_boxplot()+ facet_wrap(vars(method))+scale_fill_paletteer_d("nationalparkcolors::Acadia", name="Number of \nsampled points") +
  xlab("Number of sampled points")+ylab("Pearson correlation between suitability \nraster and prediction map of virtual species")+
  theme_minimal(base_size = 15)


p<- gridExtra::grid.arrange(model, size, points, ncol=1);rm(model, size,points)

if(!dir.exists("images/modelPerformance")) dir.create("images/modelPerformance")
ggsave(p, filename = paste0("images/modelPerformance/trueCor_byModelSizePoints.png"), dpi = 300, width = 16, height = 16)


# 8 - only PAA plots ####

for(i in c("size", "model", "points")){
  
  # Fit models dynamically using as.formula
  formula <- as.formula(paste("trueCor ~", "metric"))
  
  metrics_df <- data %>% dplyr::filter(method=="PAA")%>%
    # dplyr::filter(model %in% c("Maxent", "Lasso"))%>%
    group_by(.data[[i]]) %>%
    summarise(
      r_squared = summary(lm(formula, data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(formula, data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(formula, data = pick(everything())))[["metric"]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  p=data %>% dplyr::filter(method=="PAA") %>%
    # dplyr::filter(model %in% c("Maxent", "Lasso"))%>%
    ggplot(aes(x = .data[["metric"]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    facet_wrap(vars(!!sym(i)))+
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2), "\nRMSE = ", round(rmse, 2), "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab("Metric") +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    theme_minimal(base_size = 15)
  
  if(!dir.exists("images/PAA")) dir.create("images/PAA")
  ggsave(p, filename = paste0("images/PAA/metric_PAA_by",i,".png"), dpi = 300, width = 16, height = 8)
  
}


scale_metric <- function(value, metric) {
  
  # Definiere Baseline, Minimum, Maximum und Richtung pro Metrik
  params <- list(
    AUC   = list(baseline = 0.5, min = 0.5, max = 1, higher_better = TRUE),
    COR   = list(min = -1, max = 1, higher_better = TRUE),
    Spec  = list(min = 0, max = 1, higher_better = TRUE),
    Sens  = list(min = 0, max = 1, higher_better = TRUE),
    Kappa = list(min = -1, max = 1, higher_better = TRUE),
    PCC   = list(min = 0, max = 1, higher_better = TRUE),
    TSS   = list(baseline=0,min = -1, max = 1, higher_better = TRUE),
    PRG   = list(baseline=0.5,min = 0, max = 1, higher_better = TRUE),
    MAE   = list(min = 0, max = 1, higher_better = FALSE),
    BIAS  = list(min = -1, max = 1, higher_better = FALSE)
  )
  
  # Parameter für aktuelle Metrik holen
  if (!metric %in% names(params)) {
    stop("Metrik nicht in der Parameterliste enthalten.")
  }
  p <- params[[metric]]
  
  # Falls eine Baseline definiert ist (z.B. AUC), verwende diese als min
  if (!is.null(p$baseline)) {
    min_val <- p$baseline
  } else {
    min_val <- p$min
  }
  max_val <- p$max
  
  # Skalierung berechnen
  scaled <- (value - min_val) / (max_val - min_val)
  
  # Falls niedriger besser ist, invertieren
  if (!p$higher_better) {
    scaled <- 1 - scaled
  }
  
  # Werte auf 0-1 begrenzen
  scaled <- pmax(0, pmin(1, scaled))
  
  return(scaled)
}


