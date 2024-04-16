#' Create BUSCO plots
#'
#' This code creates BUSCO plots easier using a dataframe of summary BUSCO results
#' @param my_input_df Dataframe containing BUSCO categories C,S,D,F & M for each sample
#' @param my_total Number of total BUSCO genes in each dataset
#' @param my_title If desired provide a title
#' @param my_folder If desired provide a folder to save PDF
#' @param my_name If desired provide a name to save PDF
#' @export
my_busco_plotter <- function(my_input_df, my_total, my_title=NULL, my_folder=NULL, my_name=NULL) {

  library(reshape2)
  library(ggplot2)
  
  # my_input_df contains values for S, D, F & M categories
  # Example:
  #   Sample            C     S     D     F     M total
  # 1 Pae PA01        124   123     1     0     0   124
  # 2 Pfl UFB2        124   123     1     0     0   124
  # 3 Pfl PF5         124   124     0     0     0   124
  # 4 Pgr CPA-7       122   122     0     1     1   124
  
  ##----------------
  ## Tidy the data
  ##----------------
  colnames(my_input_df)[1] <- "my_species"
  
  my_input_df$total <- NULL
  my_input_df$C <- NULL
  my_input_df$S <- my_input_df$S/my_total*100
  my_input_df$D <- my_input_df$D/my_total*100
  my_input_df$F <- my_input_df$F/my_total*100
  my_input_df$M <- my_input_df$M/my_total*100
  
  my_input_df.m <- melt(my_input_df, id.vars = "my_species")
  colnames(my_input_df.m)[2] <- "category"
  colnames(my_input_df.m)[3] <- "my_percentage"
  
  my_input_df.m['my_values'] <- my_input_df.m$my_percentage * my_total/100
  ##----------------

  ##----------------
  ## Aesthetic Parameters
  ##----------------
  # Colors
  my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
  # Bar height ratio
  my_bar_height <- 0.75
  
  # Legend
  if (is.null(my_title)) {
    my_title <- "BUSCO Assessment Results"  
  }
  
  
  # Font
  my_family <- "sans"
  my_size_ratio <- 1
  
  my_species <- my_input_df$my_species
  
  # Code to produce the graph
  labsize = 1
  if (length(levels(my_species)) > 10){
    labsize = 0.66
  }
  ##----------------
  
  ##----------------
  ## Plot
  ##----------------
  figure <- ggplot() + 
    geom_bar(data = my_input_df.m, aes(y = my_percentage, 
                                       x = my_species, fill = category), 
             position = position_stack(reverse = TRUE), 
             stat="identity", 
             width=my_bar_height) + 
    coord_flip() + 
    theme_gray(base_size = 8) + 
    scale_y_continuous(labels = c("0","20","40","60","80","100"), 
                       breaks = c(0, 20, 40,60,80,100)) + 
    scale_fill_manual(values = my_colors,labels =c("Complete (C) and single-copy (S)  ",
                                                   "Complete (C) and duplicated (D)",
                                                   "Fragmented (F)  ",
                                                   "Missing (M)")) +   
    ggtitle(my_title) + 
    xlab("") + 
    ylab("%BUSCOs") +
    
    theme(plot.title = element_text(family=my_family, hjust=0.5, 
                                    colour = "black", 
                                    size = rel(2.2)*my_size_ratio, 
                                    face = "bold")) + 
    theme(legend.position="top",legend.title = element_blank()) + 
    theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
    theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(panel.grid.major = element_blank()) +
    theme(axis.text.y = element_text(family=my_family, 
                                     colour = "black", size = rel(1.66)*my_size_ratio)) + 
    theme(axis.text.x = element_text(family=my_family,
                                     colour = "black", size = rel(1.66)*my_size_ratio)) + 
    theme(axis.line = element_line(size=1*my_size_ratio, 
                                   colour = "black")) + 
    theme(axis.ticks.length = unit(.85, "cm")) + 
    theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
    theme(axis.ticks.x = element_line(colour="#222222")) + 
    theme(axis.ticks.length = unit(0.4, "cm")) + 
    theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
    
    guides(fill = guide_legend(override.aes = list(colour = NULL))) +
    guides(fill = guide_legend(nrow=2,byrow=TRUE))
  
  ##----------------
  
  ##----------------
  ## Add labels
  ##----------------
  for(i in rev( 1:length(levels(as.factor(my_input_df$my_species)))) ) {
    print(i)
    df_species = subset(my_input_df.m, my_species==my_input_df$my_species[i])
    print(df_species)
    total_buscos <- sum(df_species$my_values)
    figure <- figure + annotate(geom = "text", label=paste("C:", subset(df_species, category=='S')$my_values + subset(df_species, category=='D')$my_values,  
                                                           "[S:", subset(df_species, category=='S')$my_values, ", D:", subset(df_species, category=='D')$my_values," ]", 
                                                           "F:", subset(df_species, category=='F')$my_values, ", M:", subset(df_species, category=='M')$my_values, 
                                                           ", n:", total_buscos, sep=" "), 
                                y=3, x = i, 
                                size = labsize*4*my_size_ratio, 
                                colour = "black", 
                                hjust=0, family=my_family)
    
  }
  ##----------------
  
  ##----------------
  ## Return
  ##----------------
  #
  
  ## save
  if (!is.null(my_folder)) {
    pdf(file.path(my_folder, paste0(my_name, ".pdf")), paper = "A4r", width = 35, height = 12)
    print(figure)
    dev.off()
  } 
  
  return(figure)
  
}

