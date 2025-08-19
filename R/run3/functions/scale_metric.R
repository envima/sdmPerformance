
scale_metric <- function(value, metric) {
  
  # Definiere Baseline, Minimum, Maximum und Richtung pro Metrik
  params <- list(
    AUC   = list(baseline = 0.5, min = 0, max = 1, higher_better = TRUE),
    COR   = list(baseline = 0, min = -1, max = 1, higher_better = TRUE),
    trueCor   = list(baseline = 0, min = -1, max = 1, higher_better = TRUE),
    Spec  = list(min = 0, max = 1, higher_better = TRUE),
    Sens  = list(min = 0, max = 1, higher_better = TRUE),
    Kappa = list(baseline=0,min = -1, max = 1, higher_better = TRUE),
    PCC   = list(min = 0, max = 1, higher_better = TRUE),
    TSS   = list(baseline=0,min = -1, max = 1, higher_better = TRUE),
    PRG   = list(baseline=0.5,min = 0, max = 1, higher_better = TRUE),
    MAE   = list(min = 0, max = 1, higher_better = FALSE)
    
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


