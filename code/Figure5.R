#R code to generate Figure5 of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# generate graph plot -----------------------------------------------------

{
#read the saved visNetwork graph file from source_data/

conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")

# plot graph with coordinates from gephi ----------------------------------

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
## copy column "weight" to new column "value" in list "nodes"
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$weight

#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$side)
conn_graph.visn$nodes$side
#for plotting with different colors, remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

{
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)
  
visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=5, max=30),
             color = list(color='#848484', opacity=0.2),
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
    visGroups(groupname = "left_side", color=Okabe_Ito[1], shape = "dot", 
                opacity=1) %>%
  visGroups(groupname = "right_side", color=Okabe_Ito[2], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "middle", color=Okabe_Ito[6], shape = "square", 
           opacity=1) %>%
  visGroups(groupname = "fragmentum", color="white", shape = "dot", 
           opacity=0, size = 0)


visNet

#save as html
saveNetwork(visNet, "pictures/Fig5_desmo_connectome_left_right.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_left_right.html",
                    file="pictures/Fig5_desmo_connectome_left_right.png",
                    vwidth = 1500, vheight = 1500, #define the size of the browser window
                    cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}




# graph colored by segment ------------------------------------------------


#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$segment)

tail(conn_graph.visn$nodes$group)

{
#color scheme for segments
segmental_colors <- brewer.pal(6, 'Paired')
pie(rep(1,6),col=segmental_colors,segmental_colors)
  
visNet_seg <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=5, max=30),
             color = list(color='#848484', opacity=0.2),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 0.5, type = 'arrow'))) %>%
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
    visGroups(groupname = "episphere", color=segmental_colors[1], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "peristomium", color=Okabe_Ito[8], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_0", color=segmental_colors[2], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_1", color=segmental_colors[3], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_2", color=segmental_colors[4], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "segment_3", color=segmental_colors[5], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "pygidium", color=segmental_colors[6], shape = "dot", 
              opacity=1) %>%
    visGroups(groupname = "fragmentum", color="white", shape = "dot", 
              opacity=0, size = 0)
  
visNet_seg
  
#save as html
saveNetwork(visNet, "pictures/Fig5_desmo_connectome_seg.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_seg.html",
                   file="pictures/Fig5_desmo_connectome_seg.png",
                   vwidth = 1500, vheight = 1500, #define the size of the browser window
                   cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}


# color scheme for acicula ------------------------------------------------


#overwrite group value (partition) with side value or other value (for colouring)
conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$)
#need to add noto and neuropodial info to network
tail(conn_graph.visn$nodes$group)
{
visNet_seg <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=5, max=30),
           color = list(color='#848484', opacity=0.2),
           arrows = list(to = list(enabled = TRUE, 
                                   scaleFactor = 0.5, type = 'arrow'))) %>%
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
  visGroups(groupname = "episphere", color=segmental_colors[1], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "peristomium", color=Okabe_Ito[8], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "segment_0", color=segmental_colors[2], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "segment_1", color=segmental_colors[3], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "segment_2", color=segmental_colors[4], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "segment_3", color=segmental_colors[5], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "pygidium", color=segmental_colors[6], shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "fragmentum", color="white", shape = "dot", 
            opacity=0, size = 0)

visNet_acic

#save as html
saveNetwork(visNet, "pictures/Fig5_desmo_connectome_acic.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig5_desmo_connectome_acic.html",
                 file="pictures/Fig5_desmo_connectome_acic.png",
                 vwidth = 1500, vheight = 1500, #define the size of the browser window
                 cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}

}



# assemble figure --------------------------------------------------------

{


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

#read with image_read (magick package) and rotate
img_conn <- image_read("pictures/Fig1_desmo_connectome.png") %>%
image_rotate(180)

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


ggsave("figures/Figure5.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000)

}

ggsave("figures/Figure5.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000, bg = 'white')

