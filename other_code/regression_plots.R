###
ggscatter_plotRegression <- function(fitmodel, x, y, colorby, title_string="") {
  require(ggpubr)
  library(cowplot)
  p <- ggscatter(fitmodel, x, y,
                 add='reg.line',
                 conf.int = TRUE, 
                 cor.coef = TRUE,
                 palette="jco",
                 color=colorby,
                 alpha = 0.6,
                 ggtheme = theme_bw()
  ) + stat_cor(aes(color = Study), label.x=0.5) + rremove("legend")
  
  xplot <- ggboxplot(fitmodel, "Study", x, color="Study", fill="Study", palette = "jco", ggtheme = theme_bw()) + 
    rotate()
  yplot <- ggboxplot(fitmodel, "Study", y, color="Study", fill="Study", palette = "jco", ggtheme = theme_bw()) + 
    rremove("legend")
  
  #xplot <- xplot + clean_theme() + rremove("legend")
  #yplot <- yplot + clean_theme() + rremove("legend")
  
  p1 <- plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv", rel_heights = c(1,2), rel_widths = c(2,1))
  
  title_plot <- ggdraw() + draw_label(title_string, x=0.05, hjust=0) + theme(plot.margin = margin(0,0,0,7))
  
  p2 <- plot_grid(title_plot, p1, ncol = 1, rel_heights = c(0.1,1))
  
  # return plot
  return(p2)
}

## plot regression
ggplotRegression <- function (fit, Xlabel_name, Ylabel_name) {
  
  #print(fit)
  #print(Xlabel_name)
  #print(Ylabel_name)
  
  ## https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + 
    theme(plot.title = element_text(size=10)) +
    theme_classic() +
    labs(y = Ylabel_name) + labs(x = Xlabel_name)
  
}

## linear model representation
lm_plotter <- function(data, Xlabel_name, var1, var2, 
                       Ylabel_name="Counts") {
  formula_given <- paste(var1, "~", var2 )
  fit1 <- lm(formula_given, data)
  
  return(ggplotRegression(fit1, Xlabel_name, Ylabel_name))  
}