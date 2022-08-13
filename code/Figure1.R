#R code to generate panels B, C and D of Figure1 of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# read cells and connectors --------------------------------------------------------------
{
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, 
                                      fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

cirri = nlapply(read.neurons.catmaid("^anterior pair of tentacular cirri$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

muscle = nlapply(read.neurons.catmaid("^muscle$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
desmosome_connectome_non_muscle = nlapply(read.neurons.catmaid("^desmosome_connectome_non_muscle$", pid=11, 
                                                                   fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

#read muscle skids by segment with custom function (sourced from Packages_and_Connections.R)
muscle_head <- skids_by_2annotations("muscle", "episphere")
muscle_sg0 <- skids_by_2annotations("muscle", "segment_0")
muscle_sg1 <- skids_by_2annotations("muscle", "segment_1")
muscle_sg2 <- skids_by_2annotations("muscle", "segment_2")
muscle_sg3 <- skids_by_2annotations("muscle", "segment_3")
muscle_pyg <- skids_by_2annotations("muscle", "pygidium")

#read muscle skeletons and smooth
muscle_head <- nlapply(read.neurons.catmaid(muscle_head, pid=11), 
                       function(x) smooth_neuron(x, sigma=6000))
muscle_sg0 <- nlapply(read.neurons.catmaid(muscle_sg0, pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
muscle_sg1 <- nlapply(read.neurons.catmaid(muscle_sg1, pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
muscle_sg2 <- nlapply(read.neurons.catmaid(muscle_sg2, pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
muscle_sg3 <- nlapply(read.neurons.catmaid(muscle_sg3, pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
muscle_pyg <- nlapply(read.neurons.catmaid(muscle_pyg, pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))

#get the connectors
MUSconn <- connectors(muscle)

#separate synapse (prepost 1) and desmosomal (prepost 3) connectors
MUSconn1 <- MUSconn[MUSconn$prepost %in% c("1"),]
MUSconn3 <- MUSconn[MUSconn$prepost %in% c("3"),]

}

# 3d plotting -------------------------------------------------------------

{
#plot muscles coloured by segment
nopen3d() # opens a pannable 3d window
mfrow3d(1, 1)  #defines the two scenes
par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
par3d(zoom=0.52)

#panel B of Figure 2
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 

plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 

plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 

plot3d(cirri, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=1,
       add=T, alpha=0.4, col="grey") 

#color scheme for segments
segmental_colors <- brewer.pal(6, 'Paired')
pie(rep(1,6),col=segmental_colors,segmental_colors)

plot3d(muscle_head, lwd = 2, col = segmental_colors[1])
plot3d(muscle_sg0, lwd = 3, col = segmental_colors[2])
plot3d(muscle_sg1, lwd = 1, col = segmental_colors[3])
plot3d(muscle_sg2, lwd = 1, col = segmental_colors[4])
plot3d(muscle_sg3, lwd = 1, col = segmental_colors[5])
plot3d(muscle_pyg, lwd = 3, col = segmental_colors[6])

#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus.png")


#panel C of Figure 1
clear3d()
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 
plot3d(MUSconn3$x, MUSconn3$y, MUSconn3$z, add = TRUE, 
       col=hcl.colors(20000, palette='Oranges'), size=4, alpha=0.5)

#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus_desmosomes.png")

#panel D of Figure 2
clear3d()

plot3d(desmosome_connectome_non_muscle, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col=hcl.colors(length(desmosome_connectome_non_muscle), palette='Blues')) 
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 
#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus_partners.png")
close3d()
}

#clear all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# generate graph plot -----------------------------------------------------

{
#read the saved visNetwork file from source data
#this file was generated with the desmo_connectome_graph.R script
conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")
  
# plot graph with coordinates from gephi ----------------------------------

#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$class)

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
## copy column "weight" to new column "value" in list "nodes"
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$weight

#for plotting with different colors, remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()
Reds <- brewer.pal(9, 'Reds')

{
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)
#rotate the coordinates of the graph so that it is left-right symmetric
coords_rotated <- autoimage::rotate(coords, pi/2.4, pivot = c(0, 0))

visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=5, max=30),
             color = list(color='#848484', opacity=0.4),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 0.5, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(border='black'),
             opacity = 1, 
             scaling=list(min=15, max=45),
             font = list(size = 20)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',
                                       labelOnly=FALSE), 
               width = 1500, height = 1500, autoResize = FALSE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "muscle", color=Reds[8], shape = "dot", 
                opacity=1) %>%
  visGroups(groupname = "basal lamina", color="#848484", shape = "dot", 
            opacity=0.5) %>%
  visGroups(groupname = "epithelia_cell", color="#56B4E9", shape = "dot", 
           opacity=0.8) %>%
  visGroups(groupname = "circumchaetal", color="#56B4E9", shape = "dot", 
           opacity=0.8) %>%
  visGroups(groupname = "circumacicular", color="#56B4E9", shape = "dot", 
            opacity=0.8) %>%
  visGroups(groupname = "other", color="#56B4E9", shape = "dot", 
            opacity=0.8) %>%
  visGroups(groupname = "ciliated cell", color="#56B4E9", shape = "dot", 
            opacity=0.8) %>%
  visGroups(groupname = "acicula", color="#56B4E9", shape = "dot", 
            opacity=0.8) %>%
  visGroups(groupname = "chaeta", color="#56B4E9", shape = "dot", 
            opacity=0.5) %>%
  visGroups(groupname = "hemichaetal", color="#56B4E9", shape = "dot", 
            opacity=0.5) 

visNet

#save as html
saveNetwork(visNet, "pictures/Fig1_desmo_connectome.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig1_desmo_connectome.html",
                    file="pictures/Fig1_desmo_connectome.png",
                    vwidth = 1500, vheight = 1500, #define the size of the browser window
                    cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}

}

# assemble figure --------------------------------------------------------

{
panelSEM <- ggdraw() + draw_image(readPNG("images_notR/Platynereis_SEM_inverted_nolabel.png")) + 
  draw_label("head", x = 0.5, y = 0.85, size = 9) +
  draw_label("peristomium", x = 0.52, y = 0.69, size = 9) +
  draw_label("sg0", x = 0.54, y = 0.63, size = 9) +
  draw_label("sg1", x = 0.55, y = 0.54, size = 9) +
  draw_label("sg2", x = 0.56, y = 0.38, size = 9) +
  draw_label("sg3", x = 0.57, y = 0.22, size = 9) +
  draw_label("pygidium", x = 0.58, y = 0.07, size = 9) +
  draw_label("chaetae", x = 0.1, y = 0.6, size = 9) +
  draw_label("cilia", x = 0.85, y = 0.72, size = 9) +
  draw_label("Platynereis", x = 0.3, y = 0.99, size = 10, fontface = "italic") + 
  draw_label("larva", x = 0.55, y = 0.99, size = 10, fontface = "plain") +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.2, y = 0.08, size = 9)

panel_mus <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus.png")) + 
  draw_label("head", x = 0.54, y = 0.87, size = 9) +
  draw_label("sg0", x = 0.54, y = 0.69, size = 9) +
  draw_label("sg1", x = 0.52, y = 0.6, size = 9) +
  draw_label("sg2", x = 0.47, y = 0.49, size = 9) +
  draw_label("sg3", x = 0.43, y = 0.33, size = 9) +
  draw_label("pygidium", x = 0.4, y = 0.16, size = 9) +
  draw_label("cirrus sg0", x = 0.85, y = 0.72, size = 9) +
  draw_label("all muscles", x = 0.3, y = 0.99, size = 10) +
  draw_label("aciculae", x = 0.9, y = 0.28, size = 8) +
  draw_line(
    x = c(0.58, 0.8, 0.79),
    y = c(0.14, 0.28, 0.42), size = 0.3) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.82),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.82,
                   xend = 0.1,
                   yend = 0.9),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.93, size = 8) +
  draw_label("p", x = 0.1, y = 0.79, size = 8) 

panel_des <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus_desmosomes.png")) + 
  draw_label("all desmosomes", x = 0.4, y = 0.99, size = 10)

panel_des_partners <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus_partners.png")) + 
  draw_label("desmosomal partners", x = 0.45, y = 0.99, size = 10)

#read with image_read (magick package) and flip 
img_conn <- image_read("pictures/Fig1_desmo_connectome.png") %>%
  image_flop()

panel_des_conn <- ggdraw() + draw_image(img_conn) + 
  draw_label("desmosomal connectome", x = 0.45, y = 0.99, size = 10) + 
  draw_label("#desmosomes", x = 0.86, y = 0.92, size = 8) + 
  draw_label("1", x = 0.84, y = 0.88, size = 8) + 
  draw_label(max(conn_graph.visn$edges$weight), x = 0.84, y = 0.84, size = 8) +
  draw_line(x = c(0.88, 0.95), y = c(0.88, 0.88), size = 0.3, color = 'grey') +
  draw_line(x = c(0.88, 0.95), y = c(0.84, 0.84), size = 0.8, color = 'grey')

panel_desmo_EM <- ggdraw() + 
  draw_image(readPNG("images_notR/desmo_EM_Fig1.png"))

layout <- "
AAAAAABBBBBBCCCCCCDDDDDD
EEEEEEEEEEEEEEEEFFFFFFFF
"

Figure1 <- panelSEM + panel_mus + panel_des + panel_des_partners +
  panel_desmo_EM + panel_des_conn +
  plot_layout(design = layout, guides = 'collect', heights = c(1,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("figures/Figure1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000)

}

ggsave("figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000, bg = 'white')

