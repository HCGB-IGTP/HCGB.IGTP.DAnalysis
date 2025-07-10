## given a sesameQC_calcStats qcs objects, create a datframe and return it for all samples available
getDataF_sesameQC_calcStats <- function(qcs.object) {
  
  library(sesame)
  
  ## get dataframe with all info and then plot as desired
  qcs.df <- data.frame(row.names = names(qcs.object))
  for (i in 1:length(names(qcs.object))) {
    each <- sesame::sesameQC_rankStats(qc=qcs.object[[i]]) %>% as.data.frame()
    sample_ID <- names(qcs.object[i])
    each <- cbind(sample_ID, each)
    qcs.df <- rbind(qcs.df,each)
  }
  
  return(qcs.df)
}

## sesame QC step produces many quality stats, here it is a brief description
get_keys.qcs <- function() {
  keys.qcs.df <- list(
    
    "Detection" = list(
      "num_dtna" =	"N. Probes w/ Missing Raw Intensity",
      "frac_dtna" = "Probes w/ Missing Raw Intensity",
      "num_dt" = "N. Probes w/ Detection Success", 
      "frac_dt" = "% Detection Success",        
      "num_dt_mk" = "N. Detection Succ. (after masking)",
      "frac_dt_mk" = "% Detection Succ. (after masking)",  
      "num_dt_cg" = "N. Probes w/ Detection Success (cg)", 
      "frac_dt_cg" = "% Detection Success (cg)", 
      "num_dt_ch" = "N. Probes w/ Detection Success (ch)",
      "frac_dt_ch" = "% Detection Success (ch)", 
      "num_dt_rs" = "N. Probes w/ Detection Success (rs)",
      "frac_dt_rs" = "% Detection Success (rs)"),
    
    "Signal Intensity" = list(
      "mean_intensity" = "Mean sig. intensity", 
      "mean_intensity_MU" = "Mean sig. intensity (M+U)",
      "mean_ii" = "Mean sig. intensity (Inf.II)", 
      "mean_inb_grn" = "Mean sig. intens.(I.Grn IB)",  
      "mean_inb_red" = "Mean sig. intens.(I.Red IB)",  
      "mean_oob_grn" = "Mean sig. intens.(I.Grn OOB)", 
      "mean_oob_red" = "Mean sig. intens.(I.Red OOB)", 
      "na_intensity_M" = "N. NA in M (all probes)",      
      "na_intensity_U" = "N. NA in U (all probes)",
      "na_intensity_ig" = "N. NA in raw intensity (IG)",
      "na_intensity_ir" = "N. NA in raw intensity (IR)",  
      "na_intensity_ii" = "N. NA in raw intensity (II)"),
    
    "Number of Probes" = list(
      "num_probes" = "N. Probes",
      "num_probes_II" = "N. Inf.-II Probes",
      "num_probes_IR" = "N. Inf.-I (Red)",
      "num_probes_IG" = "N. Inf.-I (Grn)",
      "num_probes_cg" = "N. Probes (CG)",
      "num_probes_ch" = "N. Probes (CH)",
      "num_probes_rs" = "N. Probes (RS) "),
    
    "Color Channel" = list(
      "InfI_switch_R2R" = "N. Inf.I Probes Red -> Red",
      "InfI_switch_G2G" = "N. Inf.I Probes Grn -> Grn",
      "InfI_switch_R2G" = "N. Inf.I Probes Red -> Grn",
      "InfI_switch_G2R" = "N. Inf.I Probes Grn -> Red "),
    
    "Dye Bias" = list(
      "medR" = "Median Inf.I Intens. Red",
      "medG" = "Median Inf.I Intens. Grn",
      "topR" = "Median of Top 20 Inf.I Intens. Red",
      "topG" = "Median of Top 20 Inf.I Intens. Grn",
      "RGratio" = "Ratio of Red-to-Grn median Intens.",
      "RGdistort" = "Ratio of Top vs. Global R/G Ratios"),
    
    "Beta Value" = list(
      "mean_beta" = "Mean Beta",
      "median_beta" = "Median Beta",
      "frac_unmeth" = "% Beta < 0.3",
      "frac_meth" = "% Beta > 0.7",
      "num_na" = "N. is.na(Beta)",
      "frac_na" = "% is.na(Beta)",
      
      "mean_beta_cg" = "Mean Beta (CG)",
      "median_beta_cg" = "Median Beta (CG)",
      "frac_unmeth_cg" = "% Beta < 0.3 (CG)",
      "frac_meth_cg" = "% Beta > 0.7 (CG)",
      "num_na_cg" = "N. is.na(Beta) (CG)",
      "frac_na_cg" = "% is.na(Beta) (CG)",
      
      "mean_beta_ch" = "Mean Beta (CH)",
      "median_beta_ch" = "Median Beta (CH)",
      "frac_unmeth_ch" = "% Beta < 0.3 (CH)",
      "frac_meth_ch" = "% Beta > 0.7 (CH)",
      "num_na_ch" = "N. is.na(Beta) (CH)",
      "frac_na_ch" = "% is.na(Beta) (CH)",
      
      "mean_beta_rs" = "Mean Beta (RS)",
      "median_beta_rs" = "Median Beta (RS)",
      "frac_unmeth_rs" = "% Beta < 0.3 (RS)",
      "frac_meth_rs" = "% Beta > 0.7 (RS)",
      "num_na_rs" = "N. is.na(Beta) (RS)",
      "frac_na_rs" = "% is.na(Beta) (RS)"),
    
    ### ranking
    "rank_Intensity" = list(
      "rank_mean_intensity" = "Rank Mean sig. intensity",
      "rank_mean_ii" = "Rank Mean sig. intensity (Inf.II)", 
      "rank_mean_inb_grn" = "Rank Mean sig. intens.(I.Grn IB)", 
      "rank_mean_inb_red" = "Rank Mean sig. intens.(I.Red IB)", 
      "rank_mean_oob_grn" = "Rank Mean sig. intens.(I.Grn OOB)", 
      "rank_mean_oob_red" = "Rank Mean sig. intens.(I.Red OOB)" 
    ),
    
    "rank_NumberProbes" = list(
      "rank_num_probes" = "Rank N. Probes",  
      "rank_num_probes_II"= "Rank N. Inf.-II Probes",
      "rank_num_probes_IR" = "Rank N. Inf.-I (Red)",
      "rank_num_probes_IG" = "Rank N. Inf.-I (Grn)"  ,
      "rank_num_probes_cg"= "Rank N. Probes (CG)", 
      "rank_num_probes_ch" = "Rank N. Probes (CH)",
      "rank_num_probes_rs"= "Rank N. Probes (RS) "
    ),
    
    "rank_ColorChannel" = list(
      "rank_InfI_switch_R2R" = "Rank N. Inf.I Probes Red -> Red",
      "rank_InfI_switch_G2G" = "Rank N. Inf.I Probes Grn -> Grn",
      "rank_InfI_switch_R2G" = "Rank N. Inf.I Probes Red -> Grn", 
      "rank_InfI_switch_G2R" = "Rank N. Inf.I Probes Grn -> Red " 
    ),
    
    
    "rank_BetaValue" = list(
      "rank_mean_beta"  = "Rank Mean Beta",     
      "rank_median_beta" = "Rank Median Beta",
      "rank_frac_unmeth" = "Rank % Beta < 0.3",
      "rank_frac_meth" = "Rank % Beta > 0.7",
      "rank_num_na" = "Rank N. is.na(Beta)",        
      "rank_frac_na" = "Rank % is.na(Beta)",
      
      "rank_mean_beta_cg" = "Rank Mean Beta (CG)",
      "rank_median_beta_cg" = "Rank Median Beta (CG)",
      "rank_frac_unmeth_cg" = "Rank % Beta < 0.3 (CG)",
      "rank_frac_meth_cg" = "Rank % Beta > 0.7 (CG)",
      "rank_num_na_cg" = "Rank N. is.na(Beta) (CG)",
      "rank_frac_na_cg" = "Rank % is.na(Beta) (CG)",
      
      "rank_mean_beta_ch"  = "Rank Mean Beta (CH)",
      "rank_median_beta_ch"  = "Rank Median Beta (CH)",
      "rank_frac_unmeth_ch" = "Rank % Beta < 0.3 (CH)",
      "rank_frac_meth_ch"  = "Rank % Beta > 0.7 (CH)",
      "rank_num_na_ch" = "Rank N. is.na(Beta) (CH)", 
      "rank_frac_na_ch" = "Rank % is.na(Beta) (CH)",
      
      "rank_mean_beta_rs"= "Rank Mean Beta (RS)",
      "rank_median_beta_rs" = "Rank Median Beta (RS)",
      "rank_frac_unmeth_rs" = "Rank % Beta < 0.3 (RS)",
      "rank_frac_meth_rs" = "Rank % Beta > 0.7 (RS)",
      "rank_num_na_rs" = "Rank N. is.na(Beta) (RS)",
      "rank_frac_na_rs" = "Rank % is.na(Beta) (RS)"
    ),
    
    "rankN" = "Total Ranking"
  )
  
  return(keys.qcs.df)
}


## QC plots: check background intensity
check.background <- function(qcs.df.given, mean_oob_grn="mean_oob_grn", 
                             mean_oob_red="mean_oob_red", Sample_Name="Sample_Name", color_cat) {
  
  mean_oob_grn = ensym(mean_oob_grn)
  mean_oob_red = ensym(mean_oob_red)
  Sample_Name = ensym(Sample_Name)
  color_cat = ensym(color_cat)
  
  # plot mean_oob_grn vs. mean_oob_red
  p <- ggplot(qcs.df.given,
              aes(x = mean_oob_grn, y= mean_oob_red, label = Sample_Name, color=!!color_cat)) +
    geom_point() + geom_text(hjust = -0.1, vjust = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
    xlab('Green Background') + ylab('Red Background') +
    xlim(c(100,1000)) + ylim(c(100,1000)) + theme_light()
  
  return(p)
  
}

check.mean_intensity <- function(qcs.df.given, intensity, y_label, fill_cat){
  intensity <- ensym(intensity)
  fill_cat <- ensym(fill_cat)
  p <- ggplot(qcs.df.given) +
    geom_bar(aes(Sample_Name, !!intensity, fill=!!fill_cat), stat='identity') +
    xlab('Sample Name') + ylab(y_label) +
    ylim(0,18000) + theme_light() 
  
  print(p)
  
}

check.na_probes <- function(qcs.df.given, na2check, y_label, fill_cat){
  na2check <- ensym(na2check)
  fill_cat <- ensym(fill_cat)
  p <- ggplot(qcs.df.given) +
    geom_bar(aes(Sample_Name, !!na2check, fill=!!fill_cat), stat='identity') +
    xlab('Sample Name') + ylab(y_label) +
    theme_light() 
  
  print(p)
  
}

create_rank_plots <- function(df.given, set2plot, nameSet, pdf_file.dir) {
  
  library(fmsb)
  new.df <- df.given[,set2plot]
  
  print("Printing plots for:")
  print(nameSet)
  
  pdf_file <- file.path(pdf_file.dir, paste0(nameSet, ".pdf"))
  
  pdf(pdf_file, paper = "A4r", width = 35, height = 12)
  
  for (i in rownames(new.df)) {
    
    new.df2 <- rbind(rep(1, length(colnames(new.df[i,]))), 
                     rep(0,length(colnames(new.df[i,]))), new.df[i,])
    
    ## save in pdf
    radarchart(new.df2, axistype = 1, 
               cglcol = "black", 
               axislabcol = "black", 
               plwd = 2, 
               plty = 2, 
               pcol = rgb(0.1, 0.5, 0.5, 0.9), 
               pfcol = rgb(0.1, 0.5, 0.5, 0.5),
               cglty = 2,
               cglwd = 0.8, vlcex = 0.8) + title(paste0("Sample: ", i))
  }
  
  dev.off()
  
}

plot_identity.bar_old <- function(df.given, col2check, a.vec) {
  
  library(reshape2)
  library(ggnewscale)
  
  new.data <- df.given[,col2check]
  print(new.data)
  new.data['Sample'] <- row.names(new.data)
  print(new.data)
  
  new.data.m <- melt(new.data)
  print(new.data.m)
  
  a.vector <- create_col_palette2(columnGiven = a.vec, palette_given = "Dark2")
  
  p <- ggplot(data =new.data.m, 
              mapping = aes(x=Sample, y = value, fill=variable)) + 
    geom_bar(stat="identity") +
    theme_light() + 
    labs(color="") + 
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 8))) +
    theme(axis.text.x = element_text(color = a.vector$color.vector, face = 2))
  
  
  print(p)
}

create_col_palette2 <- function(columnGiven, palette_given = "Paired") {
  
  levels_given <- levels(as.factor(columnGiven))
  library(RColorBrewer)
  colfactors <- factor(as.factor(columnGiven), levels = levels_given)
  n_colors <- length(levels(colfactors))
  if (n_colors < 3) {
    n_colors = 3
  }
  list_colors <- brewer.pal(n = n_colors, palette_given)
  
  
  toReturn <- list(
    "color.vector" = list_colors[colfactors], 
    "named.vector" = setNames(levels_given, list_colors[1:length(levels_given)] )
  )
  
  return(toReturn)
  
  
  
  
}


plot_identity.bar <- function(df.given) {
  
  library(reshape2)
  new.data <- df.given[,c("frac_na", "frac_unmeth", "frac_meth")]
  new.data['Sample'] <- row.names(new.data)
  new.data['not_na'] <- 1 - new.data$frac_na  
  new.data['frac_middle'] <- 1 - new.data$frac_meth - new.data$frac_unmeth
  new.data['frac_meth.real'] <- new.data$not_na*new.data$frac_meth
  new.data['frac_unmeth.real'] <- new.data$not_na*new.data$frac_unmeth
  new.data['frac_middle.real'] <- new.data$not_na*new.data$frac_middle
  
  print(new.data)
  
  new.data.m <- melt(new.data[,c("Sample", "frac_meth.real", "frac_unmeth.real", "frac_middle.real", "frac_na")])
  print(new.data.m)
  
  p <- ggplot(data =new.data.m, 
              mapping = aes(x=Sample, y = value, fill=variable)) + scale_fill_brewer(palette = "Paired") +
    geom_bar(stat="identity") + 
    theme_light() + coord_flip()
  
  print(p)
}


#' Pvalue mask filtering
#' 
#' This function creates an histogram overlaid with kernel density curve
#' It uses the sesame::pOOBAH function to retrieved the pvalue detection
#' @param sample_frame SigDF generated
#' @param sample_name Sample ID
#' @param sample_dir Folder to store results
#' @param tag_name Tag to add to the files generated
pval_maskplot <- function(SigDF_given, sample_name, sample_dir, tag_name="raw_") {
  ## get results
  res_pOOBAH <- sesame::pOOBAH(SigDF_given, return.pval = TRUE)
  
  pOOBAH_df <- as.data.frame(res_pOOBAH)
  print(head(pOOBAH_df))
  
  ## Histogram overlaid with kernel density curve
  png(file.path(sample_dir, paste0(tag_name, "histogram_pvalue_filtering.png")))
  print(ggplot(pOOBAH_df, aes(x=res_pOOBAH)) + 
          geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                         binwidth=0.001,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="#FF6666") +  # Overlay with transparent density plot
          xlim(0,0.1) + geom_vline(xintercept = 0.05) + theme_light())
  dev.off()
  
  # save
  write.csv(pOOBAH_df, file = file.path(sample_dir, paste0(tag_name, 'pval_filtering_data.csv')))
}

#' Dye bias exploratory
#' 
#' Sesame Dye bias correction by matching green and red to mid point. 
#' This function compares the Type-I Red probes and Type-I Grn probes and generates 
#' and mapping to correct signal of the two channels to the middle. 
#' @param SigDF_given SigDF generated
#' @param sample_name Sample ID
#' @param sample_dir Folder to store results
#' @param tag_name Tag to add to the files generated
db_exploratory <- function(SigDF_given, sample_name, sample_dir, tag_name="raw_") {
  
  # The function takes one single SigDF and returns a SigDF with dye bias corrected.
  sset.dbNonlinear <- sesame::dyeBiasCorrTypeINorm(SigDF_given)
  
  # plot regression
  sset.dbNonlinear.plotR <- HCGB.IGTP.DAnalysis::ggscatter_plotRegression(
    data_all_given = as.data.frame(sset.dbNonlinear), x.given = 'UG', 
    y.given = 'UR', title_string = paste0("Sample: ", sample_name))
  
  png(file.path(sample_dir, paste0(tag_name, "dbNonlinear_corr.png")))
  print(sset.dbNonlinear.plotR$plot)
  dev.off()
  
  ## create qqplot
  png(file.path(sample_dir, paste0(tag_name, "dbNonlinear_qqplot.png")))
  print(qqplot(sset.dbNonlinear$UR, sset.dbNonlinear$UG,
               xlab='Type-I Red Signal', ylab='Type-I Grn Signal',
               main='Nonlinear Correction', cex=0.5))
  dev.off()
  
  # save
  save(sset.dbNonlinear, sset.dbNonlinear.plotR, 
       file = file.path(sample_dir,  paste0(tag_name, 'corr_data.RData')))
}


#######################################
