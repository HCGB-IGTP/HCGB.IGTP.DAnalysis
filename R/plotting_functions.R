#' Create PCA/UMAP plot
#'
#' This functions creates a PCA or UMAP plot given a vector of names, x and y coordinates. You can additionally
#' colour according to other vector, add title, label the points, set a different variabel to
#' include and additional parameter (size), and create circles aroung groups of row_names provided.
#' @param row_names Identifiers for each sample
#' @param x_given Variable to plot in X axis. It should be numeric and included within data_all_given
#' @param y_given Variable to plot in Y axis. It should be numeric and included within data_all_given
#' @param meta_data_col Vector for colouring
#' @param title Title to add to the plot
#' @param label_plot Add labels to the points
#' @param size TRUE/FALSE to include a numeric variable with different size points
#' @param size_var numeric variable with different size points
#' @param max.overlaps_int Number of overlaps for ggrepel
#' @param encircle Wether to include circles aroung groups of row_names provided
#' @export
plot_representation <- function(row_names, x_given, y_given, meta_data_col, title, 
                                label_plot=FALSE, size=FALSE, size_var=NULL, max.overlaps_int=50, 
                                encircle=FALSE) {
  
  
  data_plot <- data.frame(ID = row_names, 
                          x=x_given, 
                          y=y_given, 
                          col_data=meta_data_col)
  print(data_plot)
  
  if (size) {
    data_plot['size_var'] <- size_var
    data_plot <- data_plot %>% mutate(border_line=ifelse(size_var>CPOS_line, TRUE, FALSE))
    data_plot$border_line <- as.factor(data_plot$border_line)
  }
  
  pp <- ggplot(data_plot, aes(x=x, y=y))
  
  if (encircle) {
    pp <- pp + geom_encircle(aes(group=ID, col=ID), alpha=1, show.legend = TRUE, linetype=2, s_shape=0.5, expand=0)
  }
  
  if (size) {
    p <- pp + new_scale_colour() +
      geom_point(data=data_plot, aes(fill=col_data, size=size_var, colour=border_line), shape=21, show.legend = TRUE) + 
      scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white"), guide=guide_legend("Greater CPOS line"))
  } else {
    p <- pp + geom_point(aes(colour=col_data, size=21))
    
  }
  
  p =  p + ggtitle(title) + theme_classic()
  
  if (label_plot==TRUE) {
    p = p + geom_label_repel(data = data_plot, 
                             aes(x = x, y=y, label=ID), 
                             max.overlaps = max.overlaps_int)
  }
  
  print(p)
  
  return(p)
  
}

#' Shorcut to using ggscatter regression plots
#'
#' This code creates a ggscatter regression plot given two numeric variables. 
#' @param data_all_given Dataframe containing all information for all samples.
#' @param x.given Variable to plot in X axis. It should be numeric and included within data_all_given
#' @param y.given Variable to plot in Y axis. It should be numeric and included within data_all_given
#' @param title_string Title to add to the plot
#' @param x_label Title to add to the axis
#' @export
ggscatter_plotRegression <- function(data_all_given, x.given, y.given, title_string="", x_label="") {
  
  f <- as.formula(paste(x.given, y.given, sep="~"))
  print(f)
  
  fit_model <- lm(f, data_all_given)
  print (summary(fit_model))
  
  sum <- summary(fit_model)
  r.squared <- sum$r.squared
  pval <- sum$coefficients[,4][2]
  
  library(ggpubr)
  
  if (r.squared > 0.6 & pval < 0.05) {
    color4points = "red"  
  } else if (pval < 0.05) {
    color4points = "orange"  
  } else {
    color4points = "blue"  
  }
  
  
  p <- ggscatter(fit_model$model, y.given, x.given, add='reg.line', 
                 conf.int = TRUE, cor.coef = TRUE,  palette="jco", color = color4points, 
                 alpha = 0.6, ggtheme = theme_bw(), xlab = y.given) + 
    stat_cor(label.x=median(fit_model$meth))
  
  # return plot
  return(p)
}

#' Grid plots using ggscatter regression
#'
#' This code creates multiple plots using ggscatter_plotRegression for a list of variables, given
#' a coordinate (PC1, PC2, UMAP1, UMAP2...). All variables and coordinate should be included in
#' the data_all_given provided
#' @param data_all_given Dataframe containing all information for all samples.
#' @param list2test List of columns names to test. They should be valid names and should be included within data_all_given
#' @param coord.char Dimension reduction coordinate (e.g. PC1, PC2), included as a column in data_all_given
#' @export
plot.grid.ggscatterRegression <- function(list2test, data_all_given, coord.char) {
  plot_List <- list()
  for (i in list2test) {
    print(i)
    plot_List[[i]] <- ggscatter_plotRegression(data_all_given = data_all_given, 
                                               x.given = coord.char, y.given = i) 
  }
  n <- length(plot_List)
  nCol <- floor(sqrt(n))
  
  grid_plot <- arrangeGrob(grobs=plot_List, ncol = nCol)
  
  plots2return <- list(
    "plot_List" = plot_List,
    "grid_plot" = grid_plot
  ) 
  
  return(plots2return)
}

#' Shorcut to using ggscatterBoxplots
#' 
#' This code creates a ggscatter boxplot with t-test and Anova statistic by stat_compare
#' @param data_all_given Dataframe containing all information for all samples.
#' @param colname Variable name to test. It should be included within data_all_given
#' @param y.coord Dimension reduction coordinate (e.g. PC1, PC2), included as a column in data_all_given
#' @export
ggboxplot_scatter <- function(data_all_given, colName, y.coord) {
  
  library(ggpubr)
  
  comb_m <- combn(levels(factor(data_all_given[[colName]])), m = 2, simplify = TRUE)
  my_comp <- lapply(seq_len(ncol(comb_m)), function(i) comb_m[,i])
  
  gg <- ggboxplot(data_all_given, x=colName, y=y.coord) + 
    stat_compare_means(comparisons = my_comp) + 
    stat_compare_means(method = "anova", color="red") + geom_boxplot(aes(fill=colName))
  
  return(gg)
  
}

#' Grid plots using ggscatterBoxplots
#'
#' This code create multiple plots using ggboxplot_scatter for a list of variables, given
#' a coordinate (PC1, PC2, UMAP1, UMAP2...). All variables and coordinate should be included in
#' the data_all_given provided
#' @param list2test List of columns names to test. They should be valid names and should be included within data_all_given
#' @param data_all_given Dataframe containing all information for all samples.
#' @param coord.char Dimension reduction coordinate (e.g. PC1, PC2), included as a column in data_all_given
#' @export
plot.grid.ggscatterBoxplot <- function(list2test, data_all_given, coord.char) {
  plot_List <- list()
  for (i in list2test) {
    print(i)
    plot_List[[i]] <- ggboxplot_scatter(data_all_given = data_all_given, 
                                        colName = i, y.coord = coord.char) 
  }
  n <- length(plot_List)
  nCol <- floor(sqrt(n))
  
  grid_plot <- arrangeGrob(grobs=plot_List, ncol = nCol)
  
  plots2return <- list(
    "plot_List" = plot_List,
    "grid_plot" = grid_plot
  ) 
  
  return(plots2return)
}