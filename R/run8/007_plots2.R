
dataSummary$point <- factor(dataSummary$point, levels=c("40",  "80" , "120", "160", "200" ,"400"))
p1<- ggplot(dataSummary, aes(x = .data[["method"]], y = Slope)) +
  ggplot2::geom_boxplot()+xlab("")+
  theme_minimal(base_size = 24)

p2 <- ggplot(dataSummary, aes(x = .data[["method"]], y = intercept)) +
  ggplot2::geom_boxplot()+
  ylab("Intercept")+xlab("")+
  theme_minimal(base_size = 24)

p2 <- ggplot(dataSummary, aes(x = .data[["method"]], y = r2)) +
  ggplot2::geom_boxplot()+
  ylab("Intercept")+xlab("")+
  theme_minimal(base_size = 24)

dataSummary%>%dplyr::group_by(method)%>% dplyr::summarise(median = median(Slope, na.rm = TRUE))
dataSummary%>%dplyr::group_by(method)%>% dplyr::summarise(median = median(intercept, na.rm = TRUE))

p<-gridExtra::grid.arrange(p1,p2, ncol=2,nrow=1)
ggsave(p, filename = paste0("images/",nameRun, "/Intercept_Slope_boxplot.png"), dpi = 300, width = 8, height = 5)

