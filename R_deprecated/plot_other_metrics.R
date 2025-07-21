
model_PCC <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ PCC, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["PCC"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_TSS <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ TSS, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["TSS"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_Kappa <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ Kappa, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["Kappa"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

p1<- ggplot(data, aes(x = PCC, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_PCC,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")


p2<- ggplot(data, aes(x = TSS, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_TSS,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")

p3<- ggplot(data, aes(x = Kappa, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_Kappa,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")



gridExtra::grid.arrange(p1,p2,p3,nrow=3,ncol=1)#;rm(p1,p2,p3,p4,p6)




model_Spec <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ Spec, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["Spec"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_Sens <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ Sens, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["Sens"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

model_PRG <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ PRG, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["PRG"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.05,
      y = 0.95
    )
  })

p1<- ggplot(data, aes(x = Spec, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_Spec,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")


p2<- ggplot(data, aes(x = Sens, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_Sens,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")

p3<- ggplot(data, aes(x = PRG, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_PRG,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")



gridExtra::grid.arrange(p1,p2,p3,nrow=3,ncol=1)#;rm(p1,p2,p3,p4,p6)



model_MAE <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ MAE, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["MAE"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.63,
      y = 0.95
    )
  })

model_BIAS <- data %>%
  group_by(method) %>%
  do({
    fit <- lm(trueCOR ~ BIAS, data = .)
    summary_fit <- summary(fit)
    data.frame(
      slope = coef(fit)[["BIAS"]],
      r2 = summary_fit$r.squared,
      rmse = Metrics::rmse(.$trueCOR, predict(fit)),
      x = 0.63,
      y = 0.95
    )
  })


p1<- ggplot(data, aes(x = MAE, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_MAE,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")


p2<- ggplot(data, aes(x = BIAS, y = trueCOR, color = noPresencePoints)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = model_BIAS,
            aes(x = x, y = y,
                label = paste0("slope = ", round(slope, 2),
                               "\nR² = ", round(r2, 2),
                               "\nRMSE = ", round(rmse, 2))),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  ylim(0, 1) + xlim(0, 1) +
  facet_wrap(vars(method), nrow = 1) +
  scale_color_viridis(option = "D", direction = -1,name="Number of PO \npoints for testing") +
  theme_ipsum()+ylab("Pearson correlation between \nsuitability raster and prediction\n map of virtual species")



gridExtra::grid.arrange(p1,p2,nrow=3,ncol=1)#;rm(p1,p2,p3,p4,p6)

