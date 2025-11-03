#' Methods for simple_eiv objects
#'
#' Methods defined for objects returned from the error-in-variables models (e.g., dem_reg).
#'
#' @param x object of class \code{simple_eiv} from the dem_reg function.
#' @param ... further arguments passed through, see description of return value.
#'   for details.
#'   \code{\link{agree_test}}.
#' @return
#' \describe{
#'   \item{\code{print}}{Prints short summary of the error-in-variables (e.g., Deming) regression model.}
#'   \item{\code{plot}}{Returns a plot of the deming regression line, the line-of-identity, and the raw data.}
#'   \item{\code{check}}{Returns plots of the optimized residuals.}
#' }
#'
#' @name simple_eiv-methods


### methods for simple_eiv objects

#' @rdname simple_eiv-methods
#' @method print simple_eiv
#' @export
#' @source https://github.com/arcaldwell49/SimplyAgree/blob/master/R/methods.simple_agree.R (22.10.2025)


#' @rdname simple_eiv-methods
#' @method check simple_eiv
#' @export

check2 <- function(x) {
  
  df = model.frame(x$call$lm_mod)
  colnames(df) = c("y", "x")
  b0 = x$model$coef[1]
  b1 = x$model$coef[2]
  w_i = x$call$weights
  error.ratio = x$call$error.ratio
  d_i = df$y - (b0+b1*df$x)
  x_hat = df$x + (error.ratio*b1*d_i)/(1+error.ratio*b1^2)
  y_hat = df$y - (d_i)/(1+error.ratio*b1^2)
  res_x = df$x - x_hat
  res_y = df$y - y_hat
  d_sign = ifelse(d_i >= 0, 1, -1)
  opt_res = d_sign * sqrt(w_i*res_x^2 + w_i * error.ratio * res_y^2)
  avg_both = ((x_hat + y_hat )/ 2)
  
  mod <- lm(opt_res/d_sign ~ avg_both)
  
  SS <- anova(mod)$"Sum Sq"
  RegSS <- sum(SS) - SS[length(SS)]
  Chisq <- RegSS / 2
  ### Breusch-Pagan Test
  p_val_het <- pchisq(Chisq, df = 1, lower.tail = FALSE)
  df1 = data.frame(x = avg_both,
                   y = opt_res/d_sign)
  p1 = ggplot(df1,
              aes(x=x,
                  y=y)) +
    geom_point() +
    geom_smooth(se = TRUE,
                method = "loess",
                linewidth = .8,
                color ="#3aaf85",
                formula = y~x) +
    labs(y = "|Optimized Residuals|",
         x = "Average of Both Estimated Values",
         title = "Homogeneity of Residuals",
         subtitle = "Reference line should be flat and horizontal",
         caption = paste0("Heteroskedasticity", " \n",
                          "Breusch-Pagan Test: p = ",
                          signif(p_val_het,4))) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
  
  dat_norm <- na.omit(data.frame(y = opt_res))
  
  mod_res <-opt_res
  if(length(opt_res) < 5000){
    norm_test = shapiro.test(opt_res)
    norm_text = "Shapiro-Wilk Test"
  } else {
    norm_test = ks.test(opt_res, y = "pnorm",
                        alternative = "two.sided")
    norm_text = "Kolmogorov-Smirnov Test"
  }
  
 # norm_test = shapiro.test(opt_res)
  #norm_text = "Shapiro-Wilk Test"
  
  
  
  p2 = SimplyAgree:::plot_qq(x = dat_norm) +
    labs(caption = paste0("Normality", " \n",
                          norm_text, ": p = ",
                          signif(norm_test$p.value,4)))+
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
  patchwork::wrap_plots(p2, p1, ncol = 2) & patchwork::plot_annotation(
    theme = theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ))
}
