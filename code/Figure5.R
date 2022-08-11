#R code to generate Figure5 of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# generate graph plots and anatomical pictures -----------------------------------------------------

{
#read the saved visNetwork graph file from source_data/

conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")

# plot graph with coordinates from gephi ----------------------------------

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
## copy column "weight" to new column "value" in list "nodes"
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$degree

{
#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$side)
#for plotting with different colors, remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

#plot graph colored by side
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)
coords_rotated <- autoimage::rotate(coords, pi/2.4, pivot = c(0, 0))

visNet_side <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=5, max=30),
             color = list(color='#848484', opacity=0.15)) %>%
    visNodes(borderWidth=0.3, 
             color = list(border='black'),
             opacity = 1, 
             scaling=list(min=25, max=45),
             font = list(size = 0)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',
                                       labelOnly=FALSE), 
               width = 1500, height = 1500, autoResize = FALSE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "left_side", color=Okabe_Ito[1], shape = "dot", 
                opacity=1) %>%
  visGroups(groupname = "right_side", color=Okabe_Ito[2], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "middle", color=Okabe_Ito[8], size = 40, shape = "square", 
           opacity=0.8) %>%
  visGroups(groupname = "fragmentum", color="#EEEEEE", shape = "dot", 
           opacity=0, size = 0)


visNet_side

#save as html
saveNetwork(visNet_side, "pictures/Fig5_desmo_connectome_left_right.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_left_right.html",
                    file="pictures/Fig5_desmo_connectome_left_right.png",
                    vwidth = 1500, vheight = 1500, #define the size of the browser window
                    cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)

#read skeletons
left <- nlapply(read.neurons.catmaid(skids_by_2annotations("left_side", "desmosome_connectome"), pid=11), 
                function(x) smooth_neuron(x, sigma=6000))
right <- nlapply(read.neurons.catmaid(skids_by_2annotations("right_side", "desmosome_connectome"), pid=11), 
               function(x) smooth_neuron(x, sigma=6000))

middle <- nlapply(read.neurons.catmaid(skids_by_2annotations("middle", "desmosome_connectome"), pid=11), 
                 function(x) smooth_neuron(x, sigma=6000))

#plot
plot_background_ventral()
par3d(zoom=0.55)
plot3d(left, soma=TRUE, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[1]) 
plot3d(right, soma=TRUE, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[2])
plot3d(middle, soma=TRUE, lwd=3, add=T, alpha=1, col=Okabe_Ito[8])
plot3d(acicula, soma=TRUE, lwd=4, add=T, alpha=1, col="black") 
plot3d(scalebar_50um_ventral, color = "black", lw = 2)
plot3d(outline, add = TRUE, alpha = 0.017)

#make snapshot
rgl.snapshot("pictures/Fig5_left_right.png")
close3d()


}


# graph colored by segment ------------------------------------------------

#overwrite group value (partition) with side value or other value (for colouring)
{
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$segment)

#color scheme for segments
segmental_colors <- brewer.pal(6, 'Paired')
pie(rep(1,6),col=segmental_colors,segmental_colors)
  
visNet_seg <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=5, max=30),
             color = list(color='#848484', opacity=0.15)) %>%
    visNodes(borderWidth=0.3, 
             color = list(border='black'),
             opacity = 1, 
             scaling=list(min=20, max=20),
             font = list(size = 0)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',
                                       labelOnly=FALSE), 
               width = 1500, height = 1500, autoResize = FALSE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "episphere", color=segmental_colors[1], shape = "triangle", 
              opacity=1, size = 30) %>%
    visGroups(groupname = "segment_0", color=segmental_colors[2], shape = "square", 
              opacity=1) %>%
    visGroups(groupname = "segment_1", color=segmental_colors[3], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_2", color=segmental_colors[4], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_3", color=segmental_colors[5], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "pygidium", color=segmental_colors[6], shape = "square", 
              opacity=1) %>%
    visGroups(groupname = "fragmentum", color="white", shape = "dot", 
              opacity=0, size = 0)
  
visNet_seg
  
#save as html
saveNetwork(visNet_seg, "pictures/Fig5_desmo_connectome_seg.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_seg.html",
                   file="pictures/Fig5_desmo_connectome_seg.png",
                   vwidth = 1500, vheight = 1500, #define the size of the browser window
                   cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)

#read skeletons
head <- nlapply(read.neurons.catmaid(skids_by_2annotations("episphere", "desmosome_connectome"), pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
sg0 <- nlapply(read.neurons.catmaid(skids_by_2annotations("segment_0", "desmosome_connectome"), pid=11), 
                       function(x) smooth_neuron(x, sigma=6000))
sg1 <- nlapply(read.neurons.catmaid(skids_by_2annotations("segment_1", "desmosome_connectome"), pid=11), 
               function(x) smooth_neuron(x, sigma=6000))
sg2 <- nlapply(read.neurons.catmaid(skids_by_2annotations("segment_2", "desmosome_connectome"), pid=11), 
               function(x) smooth_neuron(x, sigma=6000))
sg3 <- nlapply(read.neurons.catmaid(skids_by_2annotations("segment_3", "desmosome_connectome"), pid=11), 
               function(x) smooth_neuron(x, sigma=6000))
pygidium <- nlapply(read.neurons.catmaid(skids_by_2annotations("pygidium", "desmosome_connectome"), pid=11), 
               function(x) smooth_neuron(x, sigma=6000))

#plot
plot_background_ventral()
par3d(zoom=0.55)
plot3d(head, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[1]) 
plot3d(sg0, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[2])
plot3d(sg1, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[3])
plot3d(sg2, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[4])
plot3d(sg3, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[5])
plot3d(pygidium, soma=TRUE, lwd=2, add=T, alpha=0.5, col=segmental_colors[6])
plot3d(acicula, soma=TRUE, lwd=4, add=T, alpha=1, col="black") 
plot3d(outline, add = TRUE, alpha = 0.017)

#make snapshot
rgl.snapshot("pictures/Fig5_seg.png")
close3d()


}


# color scheme for parapodia ------------------------------------------------

#overwrite group value (partition) with side value or other value (for colouring)
{
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$parapod_region)

visNet_para <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=5, max=30),
           color = list(color='#848484', opacity=0.1)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           opacity = 1, 
           scaling=list(min=30, max=55),
           font = list(size = 0)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE) %>%
  visGroups(groupname = "notopodium", color=Okabe_Ito[2], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "neuropodium", color=Okabe_Ito[1], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "non_parapodial", color="#454545", shape = "dot", 
            opacity=0.1, size = 15)

visNet_para

#save as html
saveNetwork(visNet_para, "pictures/Fig5_desmo_connectome_para.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_para.html",
                 file="pictures/Fig5_desmo_connectome_para.png",
                 vwidth = 1500, vheight = 1500, #define the size of the browser window
                 cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)

#read skeletons
notopodium <- nlapply(read.neurons.catmaid(skids_by_2annotations("notopodium", "desmosome_connectome"), pid=11), 
                  function(x) smooth_neuron(x, sigma=6000))
neuropodium <- nlapply(read.neurons.catmaid(skids_by_2annotations("neuropodium", "desmosome_connectome"), pid=11), 
                   function(x) smooth_neuron(x, sigma=6000))
acicula <- nlapply(read.neurons.catmaid(skids_by_2annotations("acicula", "desmosome_connectome"), pid=11), 
                   function(x) smooth_neuron(x, sigma=6000))

#plot
plot_background_ventral()
par3d(zoom=0.55)
plot3d(notopodium, soma=TRUE, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[2]) 
plot3d(neuropodium, soma=TRUE, lwd=2, add=T, alpha=0.4, col=Okabe_Ito[1])
plot3d(acicula, soma=TRUE, lwd=4, add=T, alpha=1, col="black") 
plot3d(outline, add = TRUE, alpha = 0.017)

#make snapshot
rgl.snapshot("pictures/Fig5_para.png")
close3d()

}

# color scheme for acicula ------------------------------------------------
#color mus, acic and circumacic
{
#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$class)

visNet_acic <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=5, max=30),
           color = list(opacity=0.5, inherit = 'both'),
           ) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='#FFFFFF'),
           opacity = 1, 
           scaling=list(min=20, max=25),
           font = list(size = 0)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE) %>%
  visGroups(groupname = "muscle", color=Reds[7], size = 25, shape = "dot", 
            opacity=0.5) %>%
  visGroups(groupname = "acicula", color=Okabe_Ito[8], size = 55, shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "circumacicular", color=Okabe_Ito[2], shape = "dot", 
            opacity=1, size = 45) %>%
  visGroups(groupname = "basal lamina", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "epithelia_cell", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "chaeta", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "other", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "circumchaetal", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "ciliated cell", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15) %>%
  visGroups(groupname = "hemichaetal", color="#FFFFFF", shape = "dot", 
            opacity=0.1, size = 15)

visNet_acic

#save as html
saveNetwork(visNet_acic, "pictures/Fig5_desmo_connectome_acic.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_acic.html",
                 file="pictures/Fig5_desmo_connectome_acic.png",
                 vwidth = 1500, vheight = 1500, #define the size of the browser window
                 cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)



#read skeletons
muscle <- nlapply(read.neurons.catmaid(skids_by_2annotations("muscle", "desmosome_connectome"), pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
circumacicular <- nlapply(read.neurons.catmaid(skids_by_2annotations("circumacicular", "desmosome_connectome"), pid=11), 
                  function(x) smooth_neuron(x, sigma=6000))
#plot
plot_background_ventral()
par3d(zoom=0.52)
plot3d(acicula, soma=TRUE, lwd=4, add=T, alpha=1, col="black") 
plot3d(circumacicular, soma=TRUE, lwd=2, add=T, alpha=0.8, col=Okabe_Ito[2]) 
plot3d(muscle, soma=TRUE, lwd=2, add=T, alpha=0.5, col=Reds[7]) 
plot3d(outline, add = TRUE, alpha = 0.017)

#make snapshot
rgl.snapshot("pictures/Fig5_mus_ac.png")
close3d()

}

}



# assemble figure --------------------------------------------------------

{

panel_left_right_3d <- ggdraw() + draw_image(readPNG("pictures/Fig5_left_right.png")) +
  draw_label("aciculae", x = 0.88, y = 0.25, size = 8) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.83, y = 0.14, size = 9) +
  draw_line(
    x = c(0.58, 0.8, 0.79),
    y = c(0.16, 0.28, 0.42), size = 0.3) +
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

panel_left_right_gr <- ggdraw() + draw_image(readPNG("pictures/Fig5_desmo_connectome_left_right.png")) + 
    draw_label("left", x = 0.3, y = 0.95, size = 8, color = Okabe_Ito[1]) + 
    draw_label("right", x = 0.7, y = 0.95, size = 8, color = Okabe_Ito[2]) + 
  draw_label("middle", x = 0.5, y = 0.95, size = 8, color = Okabe_Ito[8])

panel_seg_3d <- ggdraw() + draw_image(readPNG("pictures/Fig5_seg.png")) + 
  draw_label("head", x = 0.17, y = 0.87, size = 9, 
             color = segmental_colors[1]) + 
  draw_label("sg0", x = 0.1, y = 0.67, size = 9, color = segmental_colors[2]) +
  draw_label("sg1", x = 0.15, y = 0.58, size = 9, color = segmental_colors[3])  +
  draw_label("sg2", x = 0.13, y = 0.47, size = 9, color = segmental_colors[4]) +
  draw_label("sg3", x = 0.13, y = 0.28, size = 9, color = segmental_colors[5]) +
  draw_label("pygidium", x = 0.22, y = 0.14, size = 9, color = segmental_colors[6]) 

panel_seg_gr <- ggdraw() + draw_image(readPNG("pictures/Fig5_desmo_connectome_seg.png")) + 
  draw_label("head", x = 0.38, y = 0.90, size = 9, fontface = 'bold', 
             color = segmental_colors[1]) + 
  draw_label("sg0", x = 0.6, y = 0.9, size = 9, fontface = 'bold', color = segmental_colors[2]) +
  draw_label("sg1", x = 0.37, y = 0.79, size = 9, fontface = 'bold', color = segmental_colors[3])  +
  draw_label("sg2", x = 0.45, y = 0.55, size = 9, fontface = 'bold', color = segmental_colors[4]) +
  draw_label("sg3", x = 0.2, y = 0.2, size = 9, fontface = 'bold', color = segmental_colors[5]) +
  draw_label("pygidium", x = 0.47, y = 0.43, size = 9, fontface = 'bold', color = segmental_colors[6]) 


panel_para_3d <- ggdraw() + draw_image(readPNG("pictures/Fig5_para.png"))  

panel_para_gr <- ggdraw() + draw_image(readPNG("pictures/Fig5_desmo_connectome_para.png"))  + 
  draw_label("neuropodium", x = 0.3, y = 0.9, size = 9, fontface = 'plain', color = Okabe_Ito[1]) +
  draw_label("notopodium", x = 0.2, y = 0.72, size = 9, fontface = 'plain', color = Okabe_Ito[2])

panel_acic_3d <- ggdraw() + draw_image(readPNG("pictures/Fig5_mus_ac.png"))

panel_acic_gr <- ggdraw() + draw_image(readPNG("pictures/Fig5_desmo_connectome_acic.png"))    + 
  draw_label("aciculae", x = 0.15, y = 0.73, size = 9, color = Okabe_Ito[8]) + 
  draw_label("circumacicular", x = 0.15, y = 0.67, size = 9, color = Okabe_Ito[2]) + 
  draw_label("muscles", x = 0.15, y = 0.61, size = 9, color = Reds[7])

layout <- "
ABCD
EFGH
"

Figure5 <- panel_left_right_3d + panel_left_right_gr + panel_seg_3d + panel_seg_gr +
  panel_para_3d + panel_para_gr + panel_acic_3d+ panel_acic_gr +
  plot_layout(design = layout, guides = 'collect', widths = c(0.8, 1, 0.8, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure5.pdf", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2600, height = 1500)

}

ggsave("figures/Figure5.png", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2600, height = 1500, bg = 'white')

