##
#' Function to create single plots by control probe category
#'
#' @param data_plot Dataframe from the original ggplot in meffil qc.summary (qc.summary_given$controlmeans.summary)
#' @param vect_categories Vector containing list of categories to plot
#' @param metadata Dataframe with metadata associated to each sample
#' @param col_2_color Column to color ggplot points
#' @param col_of_interest Columns from metadata of interest
#'
#' @export
#'
plot_control_means <- function(data_plot, vect_categories, metadata, col_2_color, col_of_interest) {
  
  library(ggplot2)
  list_of_plots <- list()
  for (i_cat in vect_categories) {
    test_df <- subset(data_plot, variable==i_cat)
    test_df <- merge(metadata[,col_of_interest], test_df, by.x="Sample_Name", by.y="sample.name")
    test_df <- HCGB.IGTP.DAnalysis::df.factorizer(given.df = test_df, 
                                                  col_names.given = col_of_interest, mode = "as.factor")
    
    plot_r <- ggplot(test_df, aes_string(x="id", y="value", 
                                         colour = col_2_color)) + 
      geom_point() + stat_smooth() + ylim(c(0, max(test_df$value))) + 
      theme_light() + ggtitle(i_cat)  
    
    list_of_plots[[ i_cat ]] <- plot_r
  }
  
  return(list_of_plots)
}


#' Split vector into parts
#'
#'  Given a vector, splitted in many parts with equal number of items for later creating grid of plots
#'
#' @param vector_of_names Vector of items
#' @param max_cols Number of columns
#' @param max_rows Number of rows
#'
#' @export
#'
vector_splitter <- function(vector_of_names, max_cols = 2, max_rows=2) {
  
  # Define the number of elements in each chunk
  max_dims <- max_cols*max_rows
  max_len <- length(vector_of_names)
  
  # Split the vector into chunks and store in chunks variable
  chunks <- split(vector_of_names, ceiling(seq_along(vector_of_names) / max_dims))
  
  # create list with letters as IDs for better interpretation
  list_of_cols <- list()
  letter_id = 0
  for (i in chunks) {
    letter_id =  letter_id + 1
    letter2paste <- toupper(letters[letter_id])
    list_of_cols[[letter2paste]] <- i
  }
  return(list_of_cols)
  
}


#' Create control probe plots
#'

#' @param data_plot Dataframe from the original ggplot in meffil qc.summary (qc.summary_given$controlmeans.summary)
#' @param vect_categories Vector containing list of categories to plot
#' @param metadata Dataframe with metadata associated to each sample
#' @param col_2_color Column to color ggplot points
#' @param col_of_interest Columns from metadata of interest
#' @param tag_to_use String to include at the beggining of the PDF named generated
#' @param plots_dir Absolute path to folder to store results
#' @param max_cols Number of columns of the grid of plots
#' @param max_rows Number of rows of the grid of plots
#'
#' @export
#'
print_plot_control_means <- function(data_plot, vect_categories, col_2_color, tag_to_use,
                                     plots_dir, col_of_interest, metadata,
                                     max_cols = 2, max_rows=2){
  
  plots_control_means <- plot_control_means(data_plot = data_plot, 
                                            vect_categories = vect_categories,
                                            metadata = metadata, col_2_color = col_2_color, 
                                            col_of_interest = col_of_interest
  )
  
  cats_vector <- vector_splitter(vector_of_names = vect_categories, max_cols = max_cols, max_rows = max_rows)
  
  library(ggpubr)
  
  ## arrange plots in single ggplots
  list_of_arrange_plots <- list()
  ## save into pdf single pages
  pdf(file.path(plots_dir, paste0(tag_to_use, "control.means_by_", col_2_color,".pdf")), 
      paper = "A4r", width = 35, height = 12)
  for (i in names(cats_vector)) {
    p1 <- ggarrange(plotlist = plots_control_means[cats_vector[[i]]], 
                    ncol = max_cols, nrow = max_rows, common.legend = TRUE)
    print(p1)
    list_of_arrange_plots[[i]] <- p1
  }
  dev.off()
  
  return(plots_control_means)
  
}


#' Function to save txt information from meffil QC analysis
#'
#' Function to save as txt files all information generated during QC analysis from meffil
#' 
#' @param qc.summary_given QC summary object generated
#' @param meffil_data_dir_given Folder path to store files
#' @param tag_to_use_here String to add to each file
#'
#' @export
#'
write_qc_files <- function(qc.summary_given, meffil_data_dir_given, tag_to_use_here) {
  
  # qc.parameters
  capture.output(qc.summary_given$parameters, 
                 file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.parameters.csv")))
  
  # qc.summary$bad.cpgs
  head(qc.summary_given$bad.cpgs)
  write.csv(x = qc.summary_given$bad.cpgs, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_bad.cpgs.csv")))
  
  # qc.summary$sex.check
  head(qc.summary_given$sex.check)
  write.csv(x = qc.summary_given$sex.check, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_sex.check.csv")))
  
  ## qc.summary$bad.samples
  head(qc.summary_given$bad.samples)
  write.csv(x = qc.summary_given$bad.samples, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_bad.samples.csv")))
  
  ## qc.summary$sex.summary
  head(qc.summary_given$sex.summary$tab)
  write.csv(x = qc.summary_given$sex.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_sex.summary.csv")))
  
  ## qc.summary$meth.unmeth.summary
  head(qc.summary_given$meth.unmeth.summary$tab)
  write.csv(x = qc.summary_given$meth.unmeth.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_meth.unmeth.summary.csv")))
  
  ## qc.summary$controlmeans.summary
  head(qc.summary_given$controlmeans.summary$tab)
  write.csv(x = qc.summary_given$controlmeans.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_controlmeans.summary.csv")))
  
  ## qc.summary$sample.detectionp.summary
  head(qc.summary_given$sample.detectionp.summary$tab)
  write.csv(x = qc.summary_given$sample.detectionp.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_sample.detectionp.summary.csv")))
  
  ## qc.summary$cpg.detectionp.summary
  head(qc.summary_given$cpg.detectionp.summary$tab)
  write.csv(x = qc.summary_given$cpg.detectionp.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_cpg.detectionp.summary.csv")))
  
  ## qc.summary$sample.beadnum.summary
  head(qc.summary_given$sample.beadnum.summary$tab)
  write.csv(x = qc.summary_given$sample.beadnum.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_sample.beadnum.summary.csv")))
  
  ## qc.summary$cpg.beadnum.summary
  head(qc.summary_given$cpg.beadnum.summary$tab)
  write.csv(x = qc.summary_given$cpg.beadnum.summary$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "qc.summary_cpg.beadnum.summary.csv")))
}


#' Save QC summary plots as PDF
#'
#' @param meffil_plots_dir_given Folder path to store plots
#' @param tag_to_use_given String to add to each file
#' @param qc.summary_given QC summary object generated
#'
#' @export
#'
save_qc_plots <- function(meffil_plots_dir_given, tag_to_use_given, 
                          qc.summary_given, metadata_given, col_of_interest) {
  
  ## qc.summary$sex.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "sex.summary"), 
                                plot_given = qc.summary_given$sex.summary$graph )
  
  ## qc.summary$meth.unmeth.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "meth.unmeth.summary"), 
                                plot_given = qc.summary_given$meth.unmeth.summary$graph )
  
  ## qc.summary$controlmeans.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "controlmeans.summary"), 
                                plot_given = qc.summary_given$controlmeans.summary$graph)
  
  ###-----------------------------------------------------
  ## plot more detailed plots for control probes
  ###-----------------------------------------------------
  data_plot_controlProbes <- qc.summary_given$controlmeans.summary$graph$data
  data_plot_controlPlots <- list()
  
  data_plot_controlPlots[['study']] <- print_plot_control_means(data_plot = data_plot_controlProbes, 
                                                                vect_categories = levels(as.factor(data_plot_controlProbes$variable)),
                                                                metadata = metadata_given, col_2_color = "study", 
                                                                col_of_interest = col_of_interest, 
                                                                tag_to_use = tag_to_use_given, plots_dir = meffil_plots_dir_given, 
                                                                max_cols = 2, max_rows = 3)
  data_plot_controlPlots[['sentrix_row']] <- print_plot_control_means(data_plot = data_plot_controlProbes, 
                                                                      vect_categories = levels(as.factor(data_plot_controlProbes$variable)),
                                                                      metadata = metadata_given, col_2_color = "sentrix_row", 
                                                                      col_of_interest = col_of_interest, 
                                                                      tag_to_use = tag_to_use_given, plots_dir = meffil_plots_dir_given, 
                                                                      max_cols = 2, max_rows = 3)
  data_plot_controlPlots[['sentrix_col']] <- print_plot_control_means(data_plot = data_plot_controlProbes, 
                                                                      vect_categories = levels(as.factor(data_plot_controlProbes$variable)),
                                                                      metadata = metadata_given, col_2_color = "sentrix_col", 
                                                                      col_of_interest = col_of_interest, 
                                                                      tag_to_use = tag_to_use_given, plots_dir = meffil_plots_dir_given, 
                                                                      max_cols = 2, max_rows = 3)
  data_plot_controlPlots[['Sample_Plate']] <- print_plot_control_means(data_plot = data_plot_controlProbes, 
                                                                       vect_categories = levels(as.factor(data_plot_controlProbes$variable)),
                                                                       metadata = metadata_given, col_2_color = "Sample_Plate", 
                                                                       col_of_interest = col_of_interest, 
                                                                       tag_to_use = tag_to_use_given, plots_dir = meffil_plots_dir_given, 
                                                                       max_cols = 2, max_rows = 3)
  data_plot_controlPlots[['Slide']] <- print_plot_control_means(data_plot = data_plot_controlProbes, 
                                                                vect_categories = levels(as.factor(data_plot_controlProbes$variable)),
                                                                metadata = metadata_given, col_2_color = "Slide", 
                                                                col_of_interest = col_of_interest, 
                                                                tag_to_use = tag_to_use_given, plots_dir = meffil_plots_dir_given, 
                                                                max_cols = 2, max_rows = 3)
  ###-----------------------------------------------------
  
  ## qc.summary$sample.detectionp.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "sample.detectionp.summary"), 
                                plot_given = qc.summary_given$sample.detectionp.summary$graph)
  
  ## qc.summary$sample.beadnum.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "sample.beadnum.summary"), 
                                plot_given = qc.summary_given$sample.beadnum.summary$graph)
  
  ## qc.summary$cpg.beadnum.summary$graph
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "cpg.beadnum.summary"), 
                                plot_given = qc.summary_given$cpg.beadnum.summary$graph)
  
  ## qc.summary$cell.counts.summary$betas
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "cell.counts.summary_betas"), 
                                plot_given = qc.summary_given$cell.counts.summary$betas)
  
  ## qc.summary$cell.counts.summary$counts
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "cell.counts.summary_counts"), 
                                plot_given = qc.summary_given$cell.counts.summary$counts)
  
  ## qc.summary$genotype.summary$graphs
  HCGB.IGTP.DAnalysis::save_pdf(folder_path = meffil_plots_dir_given, 
                                name_file = paste0(tag_to_use_given, "genotype.summary"), 
                                plot_given = qc.summary_given$genotype.summary$graphs)
}

#' Save txt file for normalization objects
#'
#' @param norm.summary_given 
#' @param meffil_data_dir_given 
#' @param tag_to_use_here 
#'
#' @export
#'
write_norm_files <- function(norm.summary_given, meffil_data_dir_given, tag_to_use_here) {
  
  # qc.parameters
  capture.output(norm.summary_given$parameters, 
                 file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "parameters.csv")))
  
  # scree.plot data
  scree.plots <- names(norm.summary_given$scree.plot$tab)
  for (i in scree.plots) {
    write.csv(x = norm.summary_given$scree.plot$tab[[i]], 
              file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "scree.plots_", i, ".csv")))
  }
  
  ## control probes
  write.csv(x = norm.summary_given$control.batch$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "control.batch_tab.csv")))
  
  ## all probes
  write.csv(x = norm.summary_given$probe.batch$tab, 
            file = file.path(meffil_data_dir_given, paste0(tag_to_use_here, "probe.batch_tab.csv")))
}



#' Save and improve norm plots generated by meffil
#'
#' @param meffil_plots_dir_given 
#' @param tag_to_use_given 
#' @param norm.summary_given 
#' @param meffil_object_dir_given 
#'
#' @export
#'
save_norm_plots <- function(meffil_plots_dir_given, tag_to_use_given, norm.summary_given, meffil_object_dir_given) {
  
  ## control_probes
  control.batch.cplots <- produce_norm_plots(set_of_plots = norm.summary_given$control.batch, 
                                             name_set = "control.batch", 
                                             tag_to_use_here = tag_to_use_given, 
                                             meffil_plots_dir_given = meffil_plots_dir_given)
  
  ## all_probes
  probe.batch.cplots <- produce_norm_plots(set_of_plots = norm.summary_given$probe.batch, 
                                           name_set = "probe.batch", 
                                           tag_to_use_here = tag_to_use_given, 
                                           meffil_plots_dir_given = meffil_plots_dir_given)
  
  ## save scree.plots
  scree.plots <- norm.summary_given$scree.plot$graphs
  pdf(file.path(meffil_plots_dir_given, paste0(tag_to_use_given, "scree.plots.pdf")), 
      paper = "A4r", width = 35, height = 12)
  for (scee in names(scree.plots)) {
    print(scree.plots[[scee]])  
  }
  dev.off()
  
  ## save in directory
  save(control.batch.cplots, 
       probe.batch.cplots,
       file = file.path(meffil_object_dir_given, paste0(tag_to_use_given, "_cplots_improved.Rdata")))
  
  
}


#' Improve norm plots generated by meffil
#'
#' @param set_of_plots List of plots, either control.batch or probe.batch containing fplots, cplots and pc.plots
#' @param name_set Name of the list provided: control.batch or probe.batch
#' @param tag_to_use_here Tag to add to the file
#' @param meffil_plots_dir_given Folder to store results
#'
#' @export
#'
produce_norm_plots <- function(set_of_plots, name_set, tag_to_use_here, meffil_plots_dir_given ) {
  
  ## save fplots in pdf
  fplots_list <- set_of_plots$fplots
  pdf(file.path(meffil_plots_dir_given, paste0(tag_to_use_here, name_set, "_norm_fplots.pdf")), 
      paper = "A4r", width = 35, height = 12)
  for (i in names(fplots_list)) {
    if (get_rows(subset(fplots_list[[i]][['data']], p.value!="NA")) > 0) {
      print(fplots_list[[i]])  
    }
  }
  dev.off()
  
  ## improve cplots
  cplots_list <- set_of_plots$cplots
  cplots_list.new <- list()
  pdf(file.path(meffil_plots_dir_given, paste0(tag_to_use_here, name_set,"_cplots.pdf")), 
      paper = "A4r", width = 35, height = 12)
  for (f in names(cplots_list)) {
    data_plot <- cplots_list[[f]][['data']]
    data_plot <- data_plot %>% mutate(signif=ifelse(p.value<0.05, "True", "False"))
    data_plot$signif <- as.factor(data_plot$signif)
    p1 <- data_plot %>% 
      ggplot(mapping = aes(x=paste(x, l, sep = "."), y=estimate)) + geom_point() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, colour = signif)) +
      labs(title = paste0('Confidence Intervals for: ', f)) + 
      coord_flip() + theme_light() + xlab("Categories")
    
    cplots_list.new[[f]] <- p1
    print(p1)
  }
  dev.off()
  
  ## save pc.plots
  pc_plots <- set_of_plots$pc.plots
  pdf(file.path(meffil_plots_dir_given, paste0(tag_to_use_here, name_set, "_pc.plots.pdf")), 
      paper = "A4r", width = 35, height = 12)
  for (pcf in names(pc_plots)) {
    print(pc_plots[[pcf]])  
  }
  dev.off()
  
  
  ## return cplots improved
  return(cplots_list.new)
  
}

