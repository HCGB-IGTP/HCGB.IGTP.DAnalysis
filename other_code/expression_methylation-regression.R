##############################################
## Expression and methylation simulation
##############################################

##############################################
# Functions
##############################################

## Get parameters
get_params <- function(fit_model, comp1, comp2) {
  sum <- summary(fit_model)
  r.squared <- sum$r.squared
  pval <- sum$coefficients[,4][2]
  
  A_data <- subset(fit_model$model, Study==comp1)
  B_data <- subset(fit_model$model, Study==comp2)
  
  model_A <- lm('gene ~ meth', A_data)
  sum_A <- summary(model_A)
  r.squared_A <- sum_A$r.squared
  pval_A <- sum_A$coefficients[,4][2]
  
  model_B <- lm('gene ~ meth', B_data)
  sum_B <- summary(model_B)
  r.squared_B <- sum_B$r.squared
  pval_B <- sum_B$coefficients[,4][2]
  
  list2return <- list(
    "all_data" = list(
      "r" = cor(fit_model$model$gene, fit_model$model$meth),
      "r.squared" = r.squared,
      "pval" = pval
    ),
    "comp1_data" = list(
      "r" = cor(model_A$model$gene, model_A$model$meth),
      "r.squared" = r.squared_A,
      "pval" = pval_A
    ),
    "comp2_data" = list(
      "r" = cor(model_B$model$gene, model_B$model$meth),
      "r.squared" = r.squared_B,
      "pval" = pval_B
    )
  )
  
  return(list2return) 
}

## ggscatter regression
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

## methylation correlation
create_simulation <-function(sample_n=20, 
                              methA_min=0, methA_max=0.5,
                              exprs_meanA=1000, exprs_meanB=100,
                              slopeA=-500, slopeB=500,
                              methB_min=0.6, methB_max=1, seed=12345, t_string="Simulation") {
  
  sample_name <- rep(paste0("Sample", c(1:sample_n)))
  group1 <- factor(rep(x=c("A", "B"), each=sample_n/2))
  
  set.seed(seed)
  ## GroupA
  methA <- runif(n=sample_n/2, min = methA_min, max=methA_max)
  
  ## GroupB
  methB <- runif(n=sample_n/2, min = methB_min, max=methB_max)
  
  meth <- c(methA, methB)
  
  ## 
  epsA = rnorm(n=sample_n/2, mean = exprs_meanA, sd = sqrt(625))
  epsB = rnorm(n=sample_n/2, mean = exprs_meanB, sd = sqrt(625))
  
  print(epsA)
  print(epsB)
  
  #growth = 0 + meth*(group1=="A") + eps
  growthA = slopeA*methA + epsA
  growthB = slopeB*methB + epsB
  
  meth <- c(methA, methB)
  tmp_gene <- c(growthA, growthB)
  data_simuluation <- data.frame(sample_name, group1, meth, tmp_gene)
  data_simuluation['gene'] <- ifelse((data_simuluation$tmp_gene)<0, 0, data_simuluation$tmp_gene)
  
  rownames(data_simuluation) <- data_simuluation$sample_name
  
  print(data_simuluation)
  
  fit_beta <- lm('gene ~ meth', data_simuluation)
  print(fit_beta)
  fit_beta$model$Study <- data_simuluation$group1
  fit_beta$model$Study <- factor(fit_beta$model$Study, levels=c("A", "B"))
  bval_plot <- ggscatter_plotRegression(fitmodel = fit_beta$model, 
                                        x = "meth", y="gene", colorby = "Study", 
                                        title_string = t_string)
  print (bval_plot)
  
  return(fit_beta)
  
}
##############################################


sim1 <- create_simulation(methA_min=0, methA_max=1,
                            exprs_meanA=300, exprs_meanB=450,
                            slopeA=1, slopeB=1,
                            methB_min=0, methB_max=1, seed=123, 
                            t_string = "Simulation 1: No correlation, Diff. expression & No methylation")
get_params(sim1, "A", "B")

sim2 <- create_simulation(methA_min=0.5, methA_max=1,
                            exprs_meanA=20, exprs_meanB=20,
                            slopeA=0, slopeB=-0,
                            methB_min=0, methB_max=0.5, seed=1234,
                            t_string = "Simulation 2: No correlation, No diff expression & but diff. methylation")


sim3 <- create_simulation(methA_min=0, methA_max=0.5,
                            exprs_meanA=200, exprs_meanB=200,
                            slopeA=0, slopeB=-500,
                            methB_min=0.5, methB_max=1, seed=1234,
                            t_string = "Simulation 3: Diff. expression & partial methylation correlation (I)")


sim4 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=200, exprs_meanB=200,
                            slopeA=500, slopeB=0,
                            methB_min=0, methB_max=0.5, seed=1234, 
                            t_string = "Simulation 4: Diff. expression & partial methylation correlation (II)")

## All methylated correlation. 
sim5 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=200, exprs_meanB=200,
                            slopeA=500, slopeB=500,
                            methB_min=0, methB_max=0.5, seed=1234, 
                            t_string = "Simulation 5: High correlation, Diff. expression & high methylation correlation (I)")
sim6 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=2000, exprs_meanB=2000,
                            slopeA=-500, slopeB=-500,
                            methB_min=0, methB_max=0.5, seed=1234,
                            t_string = "Simulation 6: High correlation, Diff. expression & high methylation correlation (II)")

sim7 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=1500, exprs_meanB=2000,
                            slopeA=500, slopeB=-500,
                            methB_min=0, methB_max=0.5, seed=1234,
                            t_string = "Simulation 7: High correlation, No Diff. expression & high diff. methylation correlation (II)")

## No methylated correlation but differences
sim8 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=20, exprs_meanB=200,
                            slopeA=0, slopeB=0,
                            methB_min=0, methB_max=0.5, seed=1234,
                            t_string = "Simulation 8: No correlation, but diff. expression & methylation")

sim9 <- create_simulation2(methA_min=0.5, methA_max=1,
                            exprs_meanA=200, exprs_meanB=20,
                            slopeA=0, slopeB=0,
                            methB_min=0, methB_max=0.5, seed=1234,
                            t_string = "Simulation 9: No correlation, but diff. expression & methylation")
