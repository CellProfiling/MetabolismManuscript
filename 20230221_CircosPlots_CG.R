### author: Christian Gnann
### generates circos plots of protein multilocalization from HPA data
### used in manuscript Gnann et al: "Widespread enzyme expression heterogeneity underlies the metabolic plasticity of human cells"


library(circlize)
library(tidyverse)
library(ggsci)

date = ''
directory = "" ### where are inputfiles located
wheretosave = "" ### directory for saving of the graphs

circosplot_newcolor <- function(df, datasetstring, filepath) {
  
  colors_to_replace <- c("blue", "red", "green", "magenta", "yellow", "cyan")
  for(i in seq_along(colors_to_replace)) {
    df$color <- gsub(colors_to_replace[i],"lightgrey", df$color)
  }
  
  col = df[['color']]
  
  df$bothloc <- NULL
  df$color <- NULL
  
  loc <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  names(loc) <- c('Actin filaments','Centrosome','Cytoplasmic Aggregates', 'Cytosol','Endoplasmic reticulum',
                  'Golgi apparatus','Intermediate filaments','Microtubules','Mitochondria','Nuclear bodies',
                  'Nuclear membrane','Nuclear speckles','Nucleoli','Nucleoplasm','Plasma membrane', 'Vesicles')
  loc_order = c("Nuclear membrane","Nuclear speckles",'Nuclear bodies', "Nucleoli", "Nucleoplasm", 'Microtubules', 
                'Centrosome','Actin filaments', 'Intermediate filaments', 'Mitochondria', 'Cytosol','Cytoplasmic Aggregates',
                "Endoplasmic reticulum", "Golgi apparatus", "Plasma membrane", "Vesicles")
  ### assign colors for the compartments (nuc - blue, cyto/csk = red, and secretory = yellow)
  loc_col_comp = c('#105FAC','#105FAC','#105FAC','#105FAC','#105FAC',
              "#B9202A","#B9202A","#B9202A","#B9202A","#B9202A","#B9202A","#B9202A",
              '#1A9A2F','#1A9A2F','#1A9A2F','#1A9A2F')
  
  ### assign colors for each unique location 
  
  loc_col = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF", "#197EC0FF", "#8272B7", 
                "#46732EFF", "#71D0F5FF","#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF")
  ### make actual plot
  circos.clear()
  filename = paste(filepath, date, "_CircosPlot_",datasetstring,".pdf", sep="")
  pdf(filename) 
  #par(mar=c(6, 4, 4, 2) + 0.1)
  chordDiagram(df, order = loc_order,# annotationTrack = "grid",
               self.link = 1,
               grid.col = loc_col,
               col = col,
               transparency = 0.5,
               preAllocateTracks=list(track.height=0.5),
               link.sort = TRUE, link.largest.ontop = TRUE,)

  dev.off()
}

file_list = list.files(path = directory, all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

for (i in seq_along(file_list)) {
  filename = file_list[[i]]
  condition = str_remove(str_remove(filename, paste(date,"_CirclizeDataInput_", sep="")), ".csv")
  print(filename)
  print(condition)
  df <- read.csv(paste(directory, filename, sep="/"), header = TRUE)
  
  # run circosplot
  circosplot_newcolor(df, condition, wheretosave)
  
}
