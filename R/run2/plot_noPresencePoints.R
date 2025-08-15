library(dplyr)
library(ggplot2)

# Named list with point ranges up to 170
datasets <- list(
  "0–20 points"    = data %>% filter(noPresencePoints <= 20),
  "21–40 points"   = data %>% filter(noPresencePoints <= 40, noPresencePoints > 20),
  "41–60 points"   = data %>% filter(noPresencePoints <= 60, noPresencePoints > 40),
  "61–80 points"   = data %>% filter(noPresencePoints <= 80, noPresencePoints > 60),
  "81–100 points"  = data %>% filter(noPresencePoints <= 100, noPresencePoints > 80),
  "101–120 points" = data %>% filter(noPresencePoints <= 120, noPresencePoints > 100),
  "121–140 points" = data %>% filter(noPresencePoints <= 140, noPresencePoints > 120),
  "141–160 points" = data %>% filter(noPresencePoints <= 160, noPresencePoints > 140),
  "161–170 points" = data %>% filter(noPresencePoints <= 170, noPresencePoints > 160)
)

plot_list <- list()

for (range_label in names(datasets)) {
  
  df <- datasets[[range_label]]
  
  # Skip empty datasets to avoid errors
  if (nrow(df) == 0) next
  
  # Fit model & compute metrics
  metrics_df <- df %>%
    summarise(
      r_squared = summary(lm(as.formula(paste("trueCor ~", m)), data = cur_data()))$r.squared,
      rmse = sqrt(mean((trueCor - predict(lm(as.formula(paste("trueCor ~", m)), data = cur_data())))^2, na.rm = TRUE)),
      slope = coef(lm(as.formula(paste("trueCor ~", m)), data = cur_data()))[[m]],
      xpos = 0,
      ypos = 0.1,
      .groups = "drop"
    )
  
  # Create plot with title showing the range
  p <- ggplot(df, aes(x = .data[[m]], y = trueCor)) +
    geom_hex(bins = 140, alpha = 1) +
    scale_fill_viridis_c(name = "Point density") +
    geom_smooth(method = "lm", color = "red", fill = "#69b3a2", se = TRUE) +
    geom_text(
      data = metrics_df,
      aes(x = xpos, y = ypos, 
          label = paste0("R² = ", round(r_squared, 2),
                         "\nRMSE = ", round(rmse, 2),
                         "\nSlope = ", round(slope, 2))),
      inherit.aes = FALSE,
      size = 5,
      hjust = 0
    ) +
    xlim(0, 1) + ylim(-0.01, 1) +
    xlab(m) +
    ylab("Pearson correlation between suitability \nraster and prediction map of virtual species") +
    ggtitle(paste0("Results for ", range_label)) +
    theme_minimal(base_size = 15)
  
  # Save plot
  ggsave(filename = paste0("images/resultPlots/", m, "_", gsub(" ", "_", range_label), ".png"),
         plot = p, dpi = 300, width = 16, height = 16)
  
  plot_list[[range_label]] <- p
}



# Arrange all plots in a grid
grid_arrange <- gridExtra::grid.arrange(
  grobs = plot_list,
  ncol = 3  # 3 plots per row
)


# If you want to save the arranged grid as one image:
ggsave("images/noPresencePoints_PAA.png",
       grid_arrange,
       width = 20, height = 12, dpi = 300)
