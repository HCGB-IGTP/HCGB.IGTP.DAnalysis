#' Dataframe generator from calcStats sesameQC
#' 
#' Given a sesameQC_calcStats qcs objects, create a dataframe and return it for all samples available
#' 
#' @param qcs.object sesameQC object generated 
#'
#' @returns dataframe
#' @export
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

#' Get sesameQC keys
#'
#' sesame QC step produces many quality stats, here it is a brief description
#' @returns list
#' @export
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

#' Background intensity plotter
#'
#' Create Out-of-bag (OOB) red/green ratio plot from data frame containing values. 
#' Data frame can be generated from getDataF_sesameQC_calcStats or any other source, 
#' use appropriate name for columns conveniently 
#' 
#' @param qcs.df.given Dataframe from getDataF_sesameQC_calcStats
#' @param mean_oob_grn Column with mean oob green intensity. Default: mean_oob_grn
#' @param mean_oob_red Column with mean oob red intensity. Default: mean_oob_red
#' @param Sample_Name Column with sample ids 
#' @param color_cat Column to use for coloring samples
#' @param xlim_vec Axis limits. Default: 100-1000
#' @param ylim_vec Axis limits. Default: 100-1000
#' @returns ggplot2 object
#' @export
check.background <- function(qcs.df.given, mean_oob_grn="mean_oob_grn", 
                             mean_oob_red="mean_oob_red", Sample_Name="Sample_Name", 
                             xlim_vec = c(100,1000), ylim_vec = c(100,1000), color_cat ) {
  
  mean_oob_grn = ensym(mean_oob_grn)
  mean_oob_red = ensym(mean_oob_red)
  Sample_Name = ensym(Sample_Name)
  color_cat = ensym(color_cat)
  
  # plot mean_oob_grn vs. mean_oob_red
  p <- ggplot(qcs.df.given,
              aes(x = !!mean_oob_grn, y= !!mean_oob_red, 
                  label = !!Sample_Name, color=!!color_cat)) +
    geom_point() + 
    geom_text(hjust = -0.1, vjust = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
    xlab('Green Background') + ylab('Red Background') +
    xlim(xlim_vec) + 
    ylim(ylim_vec) + 
    theme_light()
  
  return(p)
  
}

#' Intensity plotter
#'
#' @param qcs.df.given Dataframe from getDataF_sesameQC_calcStats
#' @param intensity Intensity column name in dataframe to plot for each sample
#' @param y_label Axis Y label sentence
#' @param fill_cat Column to use for coloring samples
#' @param Sample_Name Column with sample ids 
#' 
#' @returns ggplot2 object
#' @export
check.mean_intensity <- function(qcs.df.given, intensity, y_label, fill_cat, Sample_Name="Sample_Name"){
  intensity <- ensym(intensity)
  fill_cat <- ensym(fill_cat)
  Sample_Name <- ensym(Sample_Name)
  
  p <- ggplot(qcs.df.given) +
    geom_bar(aes(!!Sample_Name, !!intensity, fill=!!fill_cat), stat='identity') +
    xlab('Sample Name') + ylab(y_label) +
    theme_light() 
  
  print(p)
}


#' Missing values NA plotter
#'
#' @param qcs.df.given Dataframe from getDataF_sesameQC_calcStats
#' @param na2check Missing values columns column name in dataframe to plot for each sample: 
#' @param y_label Axis Y label sentence
#' @param fill_cat Column to use for coloring samples
#' @param Sample_Name Column with sample ids 
#' 
#' @returns ggplot2 object
#' @export
check.na_probes <- function(qcs.df.given, na2check, y_label, fill_cat, Sample_Name="Sample_Name"){
  na2check <- ensym(na2check)
  fill_cat <- ensym(fill_cat)
  Sample_Name <- ensym(Sample_Name)
  p <- ggplot(qcs.df.given) +
    geom_bar(aes(!!Sample_Name, !!na2check, fill=!!fill_cat), stat='identity') +
    xlab('Sample Name') + ylab(y_label) +
    theme_light() 
  
  print(p)
  
}

#' Ranking plots per sample
#'
#' @param df.given Dataframe from getDataF_sesameQC_calcStats
#' @param set2plot Name of columns to use. See get_keys.qcs() for options
#' @param nameSet Name to use for saving PDF plots
#' @param pdf_file.dir Path to store files
#'
#' @export
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
    fmsb::radarchart(new.df2, axistype = 1, 
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

#' Plot identity bars
#'
#' @param df.given Dataframe from getDataF_sesameQC_calcStats
#'
#' @export
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
  
  new.data.m <- melt(new.data[,c("Sample", "frac_meth.real", 
                                 "frac_unmeth.real", "frac_middle.real", "frac_na")])
  print(new.data.m)
  
  p <- ggplot(data =new.data.m, 
              mapping = aes(x=Sample, y = value, fill=variable)) + 
    scale_fill_brewer(palette = "Paired") +
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
    data_all_given = as.data.frame(sset.dbNonlinear), 
    x.given = 'UG', 
    y.given = 'UR', 
    title_string = paste0("Sample: ", sample_name))
  
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

# https://rdrr.io/github/zwdzwd/sesame/src/R/sex.R
inferSex <- function(betas, platform=NULL) {
  hypoMALE <- c(
    "cg21983484","cg23696472","cg11673471","cg01742836","cg13574945",
    "cg08059778","cg24186901","cg26023405","cg15977272","cg13023833",
    "cg20766178","cg20455959","cg26584339","cg13130271","cg13244998",
    "cg05872808","cg21290550","cg05806018","cg07861180","cg20015269",
    "cg12576145","cg10991108","cg02333283","cg16357225","cg25206026",
    "cg20749341","cg03773146","cg04872051","cg16590821","cg09520212",
    "cg22221554","cg11152253","cg23429746","cg00813156","cg25132467",
    "cg16221895","cg09307104","cg15165114","cg18998000","cg00723973",
    "cg06041068","cg10860619","cg09514431","cg07912337","cg03334316",
    "cg17399684","cg05534333","cg23493872","cg12413138","cg05374090",
    "cg27501007","cg08855111","cg21159768","cg16488754","cg12075609",
    "cg07446674","cg01342901","cg02869694","cg12277627","cg19992190",
    "cg10717149","cg14191108","cg01869765","cg26505478","cg23685102",
    "cg02195366","cg06334238","cg02615131","cg15565409","cg15693668",
    "cg03505772","cg00845806","cg26439324","cg12935118","cg18932686",
    "cg24264679","cg08782677","cg13649400","cg06779802","cg23554546",
    "cg23951868","cg00337921","cg08479532","cg00114625","cg03391801",
    "cg22776211","cg07674503","cg22452543","cg18140045","cg15450782",
    "cg07674075","cg06510592","cg21137943","cg24479484","cg27501723",
    "cg20439892","cg18107314","cg08405463","cg09146364","cg16894263",
    "cg44822048_BC11", "cg48153389_BC11", "cg48114705_BC11", "cg48140091_BC11",
    "cg47832419_BC11", "cg47450117_BC11", "cg47728613_BC11", "cg47583295_TC11",
    "cg47476627_BC11", "cg48109634_BC11", "cg47564226_TC11", "cg47844107_BC11",
    "cg47425903_TC11", "cg47742805_BC21", "cg47855973_BC11", "cg47743423_BC11",
    "cg47906498_TC11", "cg47556267_BC11", "cg47744057_TC21", "cg48176284_BC11",
    "cg48121188_BC11", "cg48065865_BC11", "cg47748343_TC21", "cg47424030_BC21",
    "cg47744023_TC11", "cg47440985_TC11", "cg47583387_BC11", "cg47725474_TC21",
    "cg48024686_TC11", "cg47920249_BC21", "cg48114704_TC11", "cg48148849_BC21",
    "cg47742981_BC11", "cg47743136_BC11", "cg48049840_TC11", "cg48111009_TC21",
    "cg48176352_TC11", "cg47655961_BC11", "cg47856861_BC21", "cg47826283_TC11",
    "cg47901233_TC11", "cg48051845_BC11", "cg47555978_BC21", "cg47634755_BC11",
    "cg48147947_TC11", "cg47503480_BC11", "cg47740318_BC11", "cg48071477_TC11",
    "cg47643035_TC21", "cg47868567_BC11", "cg47655979_TC21", "cg47725912_BC11",
    "cg47564279_BC11", "cg48016415_TC11", "cg47656013_TC21", "cg47873187_TC21",
    "cg47438865_TC11", "cg47906673_BC11", "cg47874829_BC11", "cg47734934_BC21",
    "cg48130287_BC21", "cg47625820_BC21", "cg47505633_TC11", "cg48023062_TC21",
    "cg47744459_TC11", "cg47730002_BC11", "cg47663054_BC11", "cg47742655_BC11",
    "cg48107157_BC11", "cg48148824_TC21", "cg47634666_BC11", "cg47832434_TC11",
    "cg48057717_TC21", "cg48106464_BC11", "cg47748082_TC21", "cg47897499_BC21",
    "cg47889728_TC21", "cg47938210_TC21", "cg48176806_TC11", "cg47740347_BC11",
    "cg48021685_BC21", "cg47612856_BC11", "cg48139201_BC21", "cg48176811_BC11",
    "cg47741292_TC21", "cg47905796_TC21", "cg47643008_TC21", "cg47743984_BC11",
    "cg47795637_BC21", "cg47667056_TC11", "cg48159183_BC21", "cg48164072_TC11",
    "cg48177792_TC21", "cg47743999_TC11", "cg47471551_TC21", "cg47740813_BC21",
    "cg48157924_BC21", "cg47737568_BC11", "cg47724667_BC21", "cg47618975_BC11")
  
  hyperMALE <- c(
    "cg26359388","cg02540440","cg11049634","cg22874828","cg09182733",
    "cg01123965","cg15822015","cg05130312","cg17072671","cg22655232",
    "cg05695959","cg21010298","cg06143713","cg22759686","cg11143827",
    "cg04303560","cg11717280","cg14372935","cg05533223","cg16405492",
    "cg15765801","cg08156775","cg24183173","cg21797452","cg03161453",
    "cg10474871","cg11516614","cg18813691","cg08614574","cg08456555",
    "cg16440909","cg13326840","cg16822540","cg03801901","cg09039264",
    "cg01383599","cg14931238","cg04071644","cg22208280","cg05559023",
    "cg23317607","cg26327984","cg07801607","cg06870560","cg24156613",
    "cg04101819","cg07422795","cg14261068","cg12622895","cg09192294",
    "cg26695278","cg12653510","cg03554089","cg11166197","cg04032096",
    "cg25047306","cg07818713","cg21258987","cg07981033","cg14492530",
    "cg18157587","cg12030638","cg17498624","cg01816615","cg08723064",
    "cg05193067","cg27167763","cg15521097","cg25456959","cg16576300",
    "cg07318999","cg22417678","cg22671388","cg23644934","cg00267352",
    "cg22223709","cg23698976","cg06780606","cg13920260","cg15861835",
    "cg10039267","cg12454245","cg22067189","cg00150874","cg08401365",
    "cg13781721","cg02931660","cg01316390","cg14746118","cg21294096",
    "cg11871337","cg00408231","cg09641151","cg05226646","cg11291200",
    "cg01109660","cg23607813","cg04624564","cg07452499","cg18123612",
    "cg48211697_TC11", "cg48222828_BC21", "cg48218650_BC21", "cg48219904_BC11",
    "cg48222534_TC11", "cg48214483_TC21", "cg48222923_TC11", "cg48217358_TC12",
    "cg48217547_BC11", "cg48215035_TC11", "cg48217358_TC11", "cg48215051_TC21",
    "cg48218172_BC21", "cg48218223_TC11", "cg48296014_TC11", "cg48218620_TC11",
    "cg48213060_TC11", "cg48244014_TC11", "cg48215477_BC11", "cg48217390_BC11",
    "cg48272545_BC11", "cg48222620_BC11", "cg48309797_TC11", "cg48212920_TC11",
    "cg48218860_TC11", "cg48216374_TC11", "cg48215185_BC21", "cg48213802_BC11",
    "cg48222396_TC11", "cg48214010_BC11", "cg48222395_BC11", "cg48218465_TC21",
    "cg48215216_BC11", "cg48216938_BC11", "cg48219858_BC21", "cg48214243_BC11",
    "cg48223281_TC21", "cg48214292_BC21", "cg32022449_BC11", "cg48215159_TC21",
    "cg48222049_BC11", "cg48246403_BC21", "cg48214455_BC11", "cg48216569_TC11",
    "cg48214177_BC21", "cg48246617_BC11", "cg48301218_BC11", "cg48214011_BC11",
    "cg48215297_BC21", "cg48217555_BC21", "cg48213764_TC21", "cg48222839_TC11",
    "cg48217418_TC21", "cg48216934_BC21", "cg48250058_BC11", "cg48219493_TC21",
    "cg48222602_TC21", "cg48217485_BC12", "cg48218187_TC11", "cg48222171_TC11",
    "cg48217401_BC21", "cg48218225_BC11", "cg48222795_BC11", "cg48224019_TC11",
    "cg48217672_BC11", "cg48217626_TC21", "cg48213632_BC21", "cg48216281_TC21",
    "cg48218341_BC21", "cg48222701_BC11", "cg48218522_TC11", "cg48217489_BC11",
    "cg48212144_TC21", "cg48219215_TC21", "cg48218176_BC11", "cg48223101_BC11",
    "cg48222143_TC11", "cg48218124_BC21", "cg48218975_BC11", "cg48217449_TC21",
    "cg48222478_BC21", "cg48216323_BC21", "cg48217683_BC11", "cg48215310_TC21",
    "cg48226387_BC11", "cg48218807_BC11", "cg48213481_BC11", "cg48224372_BC11",
    "cg48217446_BC21", "cg48222402_TC11", "cg48222222_TC11", "cg48215306_BC21",
    "cg48219235_BC21", "cg48221203_TC11", "cg48216903_BC21", "cg48218631_BC21",
    "cg48220121_TC11", "cg48215553_TC11", "cg48217396_TC11", "cg48224236_BC21")
  
  platform <- sesameData_check_platform(platform, names(betas))
  if (platform != "MM285") {
    betas <- liftOver(betas, "HM450")
  }
  vals <- mean(betas[hyperMALE], na.rm = TRUE) - betas[hypoMALE]
  dd <- density(na.omit(vals))
  return(dd)
  
  if (dd$x[which.max(dd$y)] > 0.4) {
    "MALE"
  } else {
    "FEMALE"
  }
}


#######################################
