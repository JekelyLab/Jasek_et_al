#code for Video4 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely

source("code/Packages_and_Connection.R")

#read the tibble with desmosome ids, xyz coordinates and partners 1 and 2
desmo_with_partners <- read_csv("data/desmo_with_partners.csv")

celltype_names <- read.csv("data/non_neuronal_celltypes_names.csv")

#define function to read two cell classes and plot desmosomes between them and 
#skeletons from the two cell classes that are connected
plot_two_cell_types_with_desmo <- function(
  annotation1, 
  annotation2, 
  color1, 
  color2, 
  soma1, 
  soma2) {
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
  celltype1 <- nlapply(
    read.neurons.catmaid(
      unique(
        celltype1_restricted_skids
        ), 
      pid=11),
    function(x) smooth_neuron(x, sigma=6000)
    )
  celltype2 <- nlapply(
    read.neurons.catmaid(
      unique(
        celltype2_restricted_skids
        ), 
      pid=11),
    function(x) smooth_neuron(x, sigma=6000))
  
  #plot the shared desmosomes by xyz coordinates 
  plot3d(x, y, z, add = TRUE, col='red', size=6, alpha=1)
  
  #plot the two cell groups
  plot3d(
    celltype1, 
    add = TRUE, 
    soma = soma1, 
    col=color1, 
    lw=3, 
    alpha=0.4
    )
  plot3d(
    celltype2, 
    add = TRUE, 
    soma = soma2, 
    col=color2, 
    lw=3, 
    alpha=0.5
    )

}


plot_connected_desmo <- function(
  celltype1, 
  celltype2,
  celltype3, 
  textpos1, 
  textpos2, 
  textpos3) {
  #background plot
  plot_background_ventral_no_acic()
  annot1 <- celltype_names[celltype_names$Name == celltype1, "CATMAID.annotation"]
  annot2 <- celltype_names[celltype_names$Name == celltype2, "CATMAID.annotation"]
  annot3 <- celltype_names[celltype_names$Name == celltype3, "CATMAID.annotation"]
  
  plot_two_cell_types_with_desmo(
    annot1,
    annot2,
    Okabe_Ito[1], 
    Okabe_Ito[5], 
    FALSE, 
    FALSE
    )
  texts3d(textpos1, text = celltype1, col=Okabe_Ito[1], cex = 2)
  texts3d(textpos2, text = celltype2, col=Okabe_Ito[5], cex = 2)
  

  plot_two_cell_types_with_desmo(
    annot3, 
    annot2,
    Okabe_Ito[7], 
    Okabe_Ito[5],
    TRUE, 
    FALSE
    )
  texts3d(textpos3, text = celltype3, col=Okabe_Ito[7], cex = 2)
  
  #plot outline and connectome as anatomical ref
#  plot3d(connectome, add=T, alpha=0.1, lw = 1, col = "black") 
  plot3d(outline, add = T, alpha = 0.02, col = "#E2E2E2") 
  
  #snapshot with text
  for (j in 1:8) {
    rgl.snapshot(paste("videos/Video4_desmo_partners_", celltype1, celltype2, celltype3, "_first_", j, ".png", sep = ""))
  }
  
  #get the ids of text objects from the rgl window with ids3d()
  text_ids <- ids3d() %>% 
    as_tibble() %>%
    filter(type == 'text') %>%
    select(id) %>%
    pull()
  
  #remove text
  rgl.pop(type = 'shapes', id = text_ids)
  
  #full rotation with extramat - also works with many skeletons with alpha<1
  for (i in 100:200){
    nview3d("ventral", extramat=rotationMatrix((i-100)*2*pi/100, 0, 0, 1))
    #save a snapshot in the working directory
    rgl.snapshot(paste("videos/Video4_desmo_partners_", celltype1, celltype2, celltype3, "_spin_", i, ".png", sep = ""))
  }
  
  close3d()
}

# because passing a list of lists to lapply results in unexpected behavior, I have to do this manually:
plot_connected_desmo(
  "MUSlongD", 
  "MUSobP-notD", 
  "EC-circumchaetal", 
  c(49000,76000, 88000), 
  c(29000,76000, 105000), 
  c(33000,76000, 122000)
)

plot_connected_desmo(
  "circumacicular", 
  "MUSac-notA", 
  "EC", 
  c(49000,76000, 88000), 
  c(38000,76000, 70000),
  c(29000,76000, 105000)
)

plot_connected_desmo(
  "circumacicular", 
  "MUSac-notM", 
  "EC", 
  c(59000,76000, 88000), 
  c(29000,76000, 105000), 
  c(33000,76000, 122000)
)

plot_connected_desmo(
  "circumacicular", 
  "MUSac-neuDach", 
  "EC", 
  c(68000,76000, 110000), 
  c(79000,76000, 57000), 
  c(106000,76000, 150000)
)

plot_connected_desmo(
  "hemichaetal", 
  "MUSchae-neuVob", 
  "EC", 
  c(39000,76000, 52000), 
  c(32000,76000, 105000), 
  c(83000,76000, 168000)
)


#read png files and write video
library(av)
av::av_encode_video(
  paste(
    'videos/', 
    list.files(
      "videos/", 
      "Video4_desmo_partners.*.png"
      ), 
    sep = ""),
  framerate = 6,
  output = 'videos/Video4.mp4'
)


file.remove(
  paste(
    'videos/', 
    list.files(
      "videos/", 
      "Video4_.*.png"
      ), 
    sep = "")
  )
#delete individual video frames
