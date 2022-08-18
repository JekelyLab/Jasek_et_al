desmo_with_partners <- read_csv("data/desmo_with_partners.csv")

#define function to read two cell classes and plot desmosomes between them and 
#skeletons from the two cell classes that are connected
plot_two_cell_types_with_desmo <- function(annotation1, annotation2){
  annot1 <- paste("annotation:", annotation1, sep = "")
  annot2 <- paste("annotation:", annotation2, sep = "")
  celltype1_skids <- catmaid_skids(annot1, pid=11)
  celltype2_skids <- catmaid_skids(annot2, pid=11)
  
  celltype1_obP_des <- desmo_with_partners %>%
    filter(partner1 %in% celltype1_skids & partner2 %in% celltype2_skids |
             partner2 %in% celltype1_skids & partner1 %in% celltype2_skids)
  
  skids1 <- celltype1_obP_des %>%
    filter(partner1 %in% celltype1_skids) %>%
    select(partner1) %>%
    pull()
  
  skids2 <- celltype1_obP_des %>%
    filter(partner2 %in% celltype1_skids) %>%
    select(partner2) %>%
    pull()
  
  celltype1_restricted_skids <- c(skids1, skids2)
  
  skids1b <- celltype1_obP_des %>%
    filter(partner1 %in% celltype2_skids) %>%
    select(partner1) %>%
    pull()
  
  skids2b <- celltype1_obP_des %>%
    filter(partner2 %in% celltype2_skids) %>%
    select(partner2) %>%
    pull()
  
  celltype2_restricted_skids <- c(skids1b, skids2b)
  
  #get xyz coordinates
  x <- unlist(celltype1_obP_des %>% select(x))
  y <- unlist(celltype1_obP_des %>% select(y))
  z <- unlist(celltype1_obP_des %>% select(z))
  
  #read skeletons based on skids
  celltype1 <- nlapply(read.neurons.catmaid(unique(celltype1_restricted_skids), pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  celltype2 <- nlapply(read.neurons.catmaid(unique(celltype2_restricted_skids), pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  
  plot_background_ventral_no_acic()
  #plot the shared desmosomes by xyz coordinates 
  plot3d(x, y, z, add = TRUE, col=Okabe_Ito[6], size=6, alpha=1)
  
  #plot the two cell groups
  plot3d(celltype1, add = TRUE, col=Okabe_Ito[1], lw=3, alpha=0.4)
  plot3d(celltype2, add = TRUE, col=Okabe_Ito[5], lw=3, alpha=0.5)

}

plot_background_ventral_no_acic()
plot_two_cell_types_with_desmo("celltype_non_neuronal75", "celltype_non_neuronal54")

  