######################################################
#' Get p-values affected by sensitivy analysis.
#'
#' @param res_df Dataframe with results from ANCOMBC2, resulted from get_sign_data_ancombc2()
#' @param num_bins Number of bins to split the bar plot generated
#' 
#' @export
get_sensitivity_pvalues <- function(res_df, num_bins=30) {
  res_perm = res_df %>%
    dplyr::transmute(taxon = ASV,
                     p = p_val,
                     ss_pass = passed) %>%
    rowwise() %>%
    mutate(
      p_update = case_when(
        # For taxa with significant p-values that failed the sensitivity analysis,
        # we assign a random value uniformly drawn from the range [0.05, 1].
        p < 0.05 & ss_pass == FALSE ~ runif(1, min = 0.05, max = 1),
        TRUE ~ p 
      )) %>%
    ungroup()
  
  library(tidyverse)
  df_p = res_perm %>%
    dplyr::select(p, p_update) %>%
    pivot_longer(cols = p:p_update, names_to = "type", values_to = "value") 
  df_p$type = recode(df_p$type, 
                     `p` = "No Sensitivity", 
                     `p_update` = "With Sensitivity")
  
  
  fig_p <- df_p %>%
    ggplot(aes(x = value, fill = type)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = num_bins) +
    #geom_hline(yintercept = phyloseq::ntaxa(physeq.dragen)/num_bins, linetype = "dashed") +
    theme_minimal() +
    labs(title = "P-values Distribution",
         x = "P-value",
         y = "Frequency",
         fill = NULL,
         color = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom")
  fig_p
  
}

######################################################
#' Get significant results from ANCOMBC2 object
#'
#' @param res_ancombc2.obj Results from ANCOMBC2 analysis
#' @param comparison_of_interest Name of the category of interest
#' @param cat1 Name of the category used as numerator
#' @param cat2 Name of the category use as reference
#' @param q_value.given Adjusted pvalue. Default=0.05
#' @param lcf.given Absolute log fold change (LFC). Default=0.26
#'
#' @export
get_all_data_ancombc2 <- function(res_ancombc2.obj, comparison_of_interest, cat1, cat2, q_value.given=0.05, lcf.given = 0.26) {
  
  columns2get <- c("taxon",
                   "taxon",
                   paste0('lfc_', comparison_of_interest),
                   paste0('se_', comparison_of_interest),
                   paste0('W_', comparison_of_interest),
                   paste0('p_', comparison_of_interest),
                   paste0('q_', comparison_of_interest), 
                   paste0('passed_ss_', comparison_of_interest), 
                   paste0('diff_', comparison_of_interest),  
                   paste0('diff_robust_', comparison_of_interest))
  
  
  res_ancombc2.df <- res_ancombc2.obj[columns2get]
  colnames(res_ancombc2.df) <- c('taxon', 'ASV', 'lcf', 'se','W', 'p_val', 'q_value', 'passed', 'Diff_ab', 'Diff_robust')
  
  ## get data comparison
  res_ancombc2.df$Sensitivity = ifelse(res_ancombc2.df$passed, "Passed", "Not passed")
  res_ancombc2.df$association = ifelse(res_ancombc2.df$q_value < q_value.given & res_ancombc2.df$lcf > lcf.given, cat1, "Not significant")
  res_ancombc2.df$association = ifelse(res_ancombc2.df$q_value < q_value.given & res_ancombc2.df$lcf < -(lcf.given), cat2, res_ancombc2.df$association)
  
  res_ancombc2.df$ASV <- res_ancombc2.df$ASV %>% 
    stringr::str_remove_all("_Unclassified_") %>% 
    stringr::str_replace_all(pattern = "_", replacement = " ") %>% stringr::str_replace_all(pattern = "[ ]+$", replacement = "") 
  
  res_ancombc2.df 
}

######################################################
#' Create bar plot of LFC
#'
#' @param res_df Dataframe with results from ANCOMBC2, resulted from get_sign_data_ancombc2()
#' 
#' @export
bar_plot_lfc <- function(res_df) {
  
  library(scales)
  df_filter = res_df %>%
    dplyr::filter(Diff_ab == TRUE) %>% ## it is significant
    dplyr::arrange(desc(lcf)) %>%
    dplyr::mutate(direct = ifelse(lcf > 0, "Positive", "Negative"),
                  color = ifelse(Diff_robust, "black", "grey"),
                  adj = ifelse(q_value < 0.05, TRUE, FALSE))
  df_filter$ASV = factor(df_filter$ASV, levels = df_filter$ASV)
  #df_filter$direct = factor(df_filter$direct, levels = c("Positive LFC", "Negative LFC"))
  df_filter$direct = ifelse(df_filter$adj, df_filter$direct, "Not adjusted")
  
  bar_plot_lfc <- df_filter %>%
    ggplot(aes(x = ASV, y = lcf, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lcf - se, ymax = lcf + se), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log Fold change (LFC)") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1,  color = df_filter$color, face = 2)) + scale_x_discrete(labels = label_wrap(45))
  
  
  # 
  
  bar_plot_lfc
}

######################################################
#' Create Lollipop plot of LFC
#'
#' @param res_df Dataframe with results from ANCOMBC2, resulted from get_sign_data_ancombc2()
#' 
#' @export
lollipop_plot <- function(res_df) {
  
  library(scales)
  lolli_plot <- ggplot(subset(res_df, Diff_ab==TRUE), 
                       aes(y = reorder(ASV, lcf), 
                           x=lcf, color=association, 
                           shape=Sensitivity)) + 
    geom_point(aes(size=q_value)) + 
    scale_size(trans = 'reverse') + 
    theme_bw() +   
    scale_shape_manual(values = c(1,16)) + 
    geom_segment(aes(
      xend = 0, 
      y = reorder(ASV, lcf),
      yend = reorder(ASV, lcf)), color = "darkgrey") +
    geom_vline(xintercept = 0, size=0.3, linetype = 'dashed', alpha = 0.25) + 
    xlab("LFC") +ylab(NULL) + 
    labs(title = "Lollipop plot") +
    scale_y_discrete(labels = label_wrap(45))
  
  ## return
  lolli_plot
  
}

######################################################
#' Create volcano plot of LFC vs. q_value
#'
#' @param res_df Dataframe with results from ANCOMBC2, resulted from get_sign_data_ancombc2()
#' 
#' @export
volcano_ancombc2 <- function(res_df) {
  vplot <- ggplot(res_df, aes(x = as.numeric(lcf), 
                              y = -log10(as.numeric(q_value)), 
                              color = association, 
                              shape=Sensitivity,
                              label=ASV)) +
    geom_point() + labs(x = "LFC", 
                        y = "-log10(q-value)", 
                        title = "Volcano Plot") + 
    scale_shape_manual(values = c(1,16)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + ggrepel::geom_text_repel() + theme_bw() 
  vplot 
}

######################################################
#' Create all DE plots 
#'
#' @param res_df Dataframe with results from ANCOMBC2, resulted from get_sign_data_ancombc2()
#' 
#' @export
DE_plots_ancombc2 <- function(res_df) {
  
  ## volcano plot
  library(ggplot2)
  library(patchwork)
  
  volcano <- volcano_ancombc2(res_df) 
  lollipot  <-  lollipop_plot(res_df)
  barplot.here  <- bar_plot_lfc(res_df = res_df)
  sensitivity  <- get_sensitivity_pvalues(res_df = res_df)
  
  arranged_plot <- (volcano / sensitivity) | lollipot
  
  arranged_plot <- arranged_plot + plot_annotation(tag_levels = 'A')
  
  return(list("volcano"=volcano, 
              "lollipot"=lollipot,
              "barplot"=barplot.here, 
              "sensitivity"=sensitivity,
              "arranged_plot" = arranged_plot)
  )
  
}

######################################################
#' Create string with taxonomy information
#'
#' @param tax2use ID resulted from ANCOMBC2 analysis
#' @param tax_table2use Original taxonomic information dataframe
#' @param retString Return string or dataframe
#' @param tax_level Taxonomic level to include, as used in ANCOMBC2 analysis
#'
#' @export
get_string_taxonomy <- function(tax2use, tax_table2use, retString=TRUE, tax_level="Genus") {
  
  ## Get the levels to include
  tax_levels=c("Domain","Phylum","Class",
               "Order","Family","Genus", "Species")
  
  tax_levels <- tax_levels[1:match(tax_level, tax_levels)]
  
  string2return <- ""
  if (stringr::str_detect(string = tax2use, pattern = "Unclassified")) {
    if (retString) {
      return(tax2use)
    } else {
      df2use = data.frame(row.names = tax_levels, values=c(rep("", 6)))
      df2use[1,'values'] <- tax2use
      df2use <- t(df2use) %>% data.frame()
      return(df2use)
    }
  }
  
  if (tax2use %in% tax_table2use[[tax_level]]) {
    df2use <- tax_table2use[tax_table2use[[tax_level]]==tax2use,]
    df2use <- df2use[,tax_levels]
    rownames(df2use) <- 1:nrow(df2use)
    df2use <- unique(df2use)
    
  } else {
    
    strin2use <- tax2use %>% stringr::str_replace_all("[_]+", " -") %>% 
      stringr::str_split("-") %>% unlist()
    
    if (length(strin2use) == length(tax_levels)) {
      df2use <- data.frame(strin2use, row.names = tax_levels)
    } else {
      diff_len <- length(tax_levels) - length(strin2use)
      
      if (!diff_len>0) {
        if (retString) {
          return(tax2use)
        } else {
          df2use = data.frame(row.names = tax_levels, values=c(rep("", 6)))
          df2use[1,'values'] <- tax2use
          df2use <- t(df2use) %>% data.frame()
          return(df2use)
        }
      }
      
      strin2use <- c(strin2use, rep(x = "", diff_len))
      df2use <- data.frame(strin2use, row.names = tax_levels)
      
    }
    df2use <- t(df2use) %>% as.data.frame()
    df2use[df2use==""]="Unclassified"
  }
  
  string2return = paste(paste(col(df2use, TRUE), as.matrix(df2use), sep = ": "), collapse = "; ")
  
  if (retString) {
    return(string2return)  
  } else {
    return(df2use)
  }
  
}

######################################################
#' Get significant summary results for ANCOMBC2
#'
#' @param ancombc_res ANCOMBC2 object results
#' @param filter_passed TRUE/FALSE to include filter passed results
#' @param filter_diff TRUE/FALSE to include filter different results
#' @param min_abs_lfc Log Fold Change cutoff to include
#' 
#' @export
signif_ancombc <- function(ancombc_res, filter_passed = TRUE, filter_diff = TRUE, min_abs_lfc = 0) {

  ## Original idea: https://github.com/adrientaudiere/MiscMetabar/blob/4cd6ffc736d6dd0ab166a8147dab157d44d77b89/R/beta_div_test.R#L848

    
  signif_ancombc_res <- ancombc_res
  clnames <- colnames(signif_ancombc_res)
  
  name_modality <- gsub("passed_ss", "", 
                        clnames[grepl("passed_ss", clnames) & !grepl("Intercept", clnames)])
  
  print("..........................")
  print("name_modality")
  print(name_modality)
  print(paste0("passed_ss", name_modality))
  print("..........................")
  
  ## all filtered
  if (filter_passed) {
    signif_ancombc_res <- signif_ancombc_res %>%
      filter(Reduce(`|`, select(., c(paste0("passed_ss", name_modality)))))
  }
  
  if (filter_diff) {
    signif_ancombc_res <- signif_ancombc_res %>% 
      filter(Reduce(`|`, select(., c(paste0("diff", name_modality)))))
  }
  
  signif_ancombc_res.pos <- signif_ancombc_res %>%
    filter(Reduce(`|`, select(., c(paste0("lfc", name_modality)))) > min_abs_lfc)
  
  signif_ancombc_res.neg <- signif_ancombc_res %>%
    filter(Reduce(`|`, select(., c(paste0("lfc", name_modality)))) < -(min_abs_lfc))
  
  res_list <- list()
  res_list[['all']] = list("pos"=signif_ancombc_res.pos,
                           "neg"=signif_ancombc_res.neg)
  
  for (n in name_modality) {
    signif_ancombc_res <- ancombc_res
    if (filter_passed) {
      signif_ancombc_res <- signif_ancombc_res %>%
        filter(Reduce(`|`, select(., paste0("passed_ss", n))))
    }
    
    if (filter_diff) {
      signif_ancombc_res <- signif_ancombc_res %>% 
        filter(Reduce(`|`, select(., paste0("diff", n))))
    }
    
    ## get positive
    signif_ancombc_res.pos <- signif_ancombc_res %>%
      filter(Reduce(`|`, select(., paste0("lfc", n))) > min_abs_lfc)
    
    
    signif_ancombc_res.neg <- signif_ancombc_res %>%
      filter(Reduce(`|`, select(., paste0("lfc", n))) < -(min_abs_lfc))
    
    res_list[[n]] =  list("pos"=signif_ancombc_res.pos,
                          "neg"=signif_ancombc_res.neg)
    
  }
  
  return(res_list)
}

######################################################
#' Get results for ANCOMBC2 comparison of interest
#'
#' @param path_given Main folder to store results
#' @param comp.given_name Name of the comparison of interest, to use during saving results
#' @param comp.interest Name of the comparison of interest as reported in the ANCOMBC2 results
#' @param res_ancombc2_results.df Res Dataframe returned by ANCOMBC2 or by signif_ancombc()
#' @param res_list List of results returned by res_ancombc2_RefC()
#'
#' @export
get_results_Ancombc <- function(path_given, comp.given_name, comp.interest, res_ancombc2_results.df) {
  
  dir_path.name <- file.path(path_given, comp.given_name)
  dir.create(dir_path.name)
  dir(dir_path.name)
  
  results_comp <- HCGB.IGTP.DAnalysis::get_comparison_resultsNames(paste0('test_', comp.given_name))
  
  all_res_df <- get_all_data_ancombc2(res_ancombc2.obj = res_ancombc2_results.df, 
                                                comparison_of_interest = comp.interest, 
                                                cat1 = results_comp$cmp1, cat2 = results_comp$cmp2)
  
  DE_plots.list <- DE_plots_ancombc2(all_res_df)
  
  ## Save in PDF
  HCGB.IGTP.DAnalysis::save_multi_pdf(folder_path = dir_path.name, 
                                      name_file = paste0(comp.given_name, "-DEplots"), list_plots = DE_plots.list)
  
  sign_res_df <- subset(all_res_df, association != "Not significant")
  results_list2save <- list("all" = all_res_df,
                            "sign_results"= sign_res_df,
                            "robust" = subset(all_res_df, Diff_robust==TRUE))
  
  HCGB.IGTP.DAnalysis::save_woorkbook_openxlsx(file_path = file.path(dir_path.name, paste0(comp.given_name, "-results.xlsx")), 
                                               list_df = results_list2save)
  
  write.csv(all_res_df, file = file.path(dir_path.name, paste0(comp.given_name, "-all_results.csv")))
  write.csv(sign_res_df, file = file.path(dir_path.name, paste0(comp.given_name, "-sign_results.csv")))
  
  save(DE_plots.list, results_list2save,
       file = file.path(dir_path.name, paste0(comp.given_name, ".RData")))
  
  return(list("DE_plots.list" = DE_plots.list, 
              "results_list2save" = results_list2save))
  
}

