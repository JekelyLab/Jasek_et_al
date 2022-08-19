#code for Figure6 figure supplement5 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely

source("code/Packages_and_Connection.R")

#read the tibble with desmosome ids, xyz coordinates and partners 1 and 2
desmo_with_partners <- read_csv("data/desmo_with_partners.csv")

#read skids of connectome neurons to plot as anatomical ref
connectome_skids <- catmaid_skids("^connectome_neuron$", pid=11)
#read a random subset of connectome skeletons for speed up
connectome <- nlapply(read.neurons.catmaid(sample(connectome_skids, 300), pid=11),
                      function(x) smooth_neuron(x, sigma = 6000))

#define function to read two cell classes and plot desmosomes between them and 
#skeletons from the two cell classes that are connected
plot_two_cell_types_with_desmo <- function(annotation1, annotation2, 
                                           color1, color2, 
                                           soma1, soma2){
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
  
  #plot the shared desmosomes by xyz coordinates 
  plot3d(x, y, z, add = TRUE, col='red', size=6, alpha=1)
  
  #plot the two cell groups
  plot3d(celltype1, add = TRUE, soma = soma1, col=color1, lw=3, alpha=0.4)
  plot3d(celltype2, add = TRUE, soma = soma2, col=color2, lw=3, alpha=0.5)

}

plot_background_ventral_no_acic()
#plot MUSlongD and MUSobP_notD
plot_two_cell_types_with_desmo("celltype_non_neuronal75", 
                               "celltype_non_neuronal54",
                               Okabe_Ito[1], Okabe_Ito[5], 
                               FALSE, FALSE)
texts3d(49000,76000, 88000, text = "MUSlongD", col='black', cex = 2)
texts3d(29000,76000, 105000, text = "MUSobP-notD", col='black', cex = 2)
texts3d(33000,76000, 122000, text = "EC-circumchaetal", col='black', cex = 2)

#plot MUSobP_notD and EC-circumchaetal
plot_two_cell_types_with_desmo("celltype_non_neuronal28", 
                               "celltype_non_neuronal54",
                               Okabe_Ito[7], Okabe_Ito[5],
                               TRUE, FALSE)

#plot outline and connectome as anatomical ref
plot3d(connectome, add=T, alpha=0.1, lw = 1, col = "black") 
plot3d(outline, add = T, alpha = 0.02, col = "#E2E2E2") 

#snapshot with text
rgl.snapshot("pictures/Fig6_fig_suppl5_1.png")

#do the same for other muscle partners


