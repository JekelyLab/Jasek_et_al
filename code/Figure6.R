#code for Figure6 of the Jasek et al 2021 desmosomal connectome paper
#quantification and visualisation of node degrees in the desmosomal connectome 
#Gaspar Jekely

source("code/Packages_and_Connection.R")

#read the saved igraph format graph file from supplements/
desmo_conn_graph <- readRDS("source_data/Figure3_source_data1.rds")
V(desmo_conn_graph)$weighted_degree <- V(desmo_conn_graph)$degree

V(desmo_conn_graph)$degree <- degree(
  desmo_conn_graph,
  v = V(desmo_conn_graph),
  mode = "all",
  loops = TRUE,
  normalized = FALSE
)

node_name_and_degree_high_only <- as_tibble(desmo_conn_graph %>%  
  as_tbl_graph() %>%
  filter(weighted_degree > 10) %>%
  select(name, weighted_degree))

# 3d plotting  ------------------------------------------------------------

{
plot_background_ventral()
par3d(zoom=0.55)

#plot skeletons with transparency inversely proportional to weighted degree
for (i in 1:length(node_name_and_degree_high_only$name)) {
  skeleton <- nlapply(read.neurons.catmaid(node_name_and_degree_high_only$name[i], pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
  plot3d(skeleton, soma = TRUE, lwd = 1, add=T, 
         alpha=node_name_and_degree_high_only$weighted_degree[i]/max(node_name_and_degree_high_only$weighted_degree),
         col=Reds[8])
}

plot3d(outline, add = TRUE, alpha = 0.017)
plot3d(acicula, soma=TRUE, lwd=2, add=T, alpha=0.7, col="black") 
plot3d(scalebar_50um_ventral, color = "black", lw = 2)

#make snapshot
rgl.snapshot("pictures/Fig6_desmo_connectome_highest_weight_ventral.png")

#remove scalebar for lateral view
rgl.pop(type = "shapes")

#lateral view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.55)
rgl.snapshot("pictures/Fig6_desmo_connectome_highest_weight_right.png")
close3d()
}

# plot graph with nodes coloured by weighted degree ------------------------

{
#read the saved visNetwork graph file from source_data/
conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
## copy column "weight" to new column "value" in list "nodes"
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$weighted_degree

#remove colour and group info
conn_graph.visn$nodes$color <- c()
conn_graph.visn$nodes$group <- c()
#define node labels
conn_graph.visn$nodes$label <- conn_graph.visn$nodes$CATMAID_name
#assign colour on a red to white gradient depending on node degree
colors <- c('white', Reds)

conn_graph.visn$nodes$color <- colors[round(conn_graph.visn$nodes$weighted_degree/max(conn_graph.visn$nodes$weighted_degree)*9)+1]

#coordinates matrix
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)
#rotate the coordinates of the graph so that it is left-right symmetric
coords_rotated <- autoimage::rotate(coords, pi/2.4, pivot = c(0, 0))

#filter labels, only keep those with high weighted degree
conn_graph.visn$nodes$label <- tibble(conn_graph.visn$nodes) %>%
  select(label, weighted_degree) %>%
  mutate(label = ifelse(weighted_degree > 45, label, '')) %>%
  select(label) %>%
  pull

visNet_degree <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=5, max=30),
           color = list(source='both', opacity=0.2)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           opacity = 1, 
           scaling=list(min=1, max=55),
           font = list(size = 70)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE)

visNet_degree

#save as html
saveNetwork(visNet_degree, "pictures/Fig6_desmo_connectome_degree.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig6_desmo_connectome_degree.html",
                 file="pictures/Fig6_desmo_connectome_degree.png",
                 vwidth = 1500, vheight = 1500, #define the size of the browser window
                 cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)

}

# make plot of node weighted degree and degree  ----------------------------

{
degree_plot <- as_tibble(desmo_conn_graph %>%  
  as_tbl_graph()) %>%
  filter(class != "other") %>%
  filter(class != "basal lamina") %>%
  select(class, degree) %>%
  ggplot(aes(x = class, y = degree)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 0.6, aes(color = degree, alpha = degree/84)) +
  scale_colour_gradient(low = Reds[1], high = Reds[9], 
                        guide = "colorbar", 
                        guide_legend('degree')) +
  labs(y = 'degree', x = 'cell class') +
  theme_half_open() + 
  guides(alpha = "none", color = "colorbar") +
  theme(panel.grid.major = element_blank(),
                      panel.border = element_blank(),
                      axis.text.x = element_text (angle = 90, hjust = 1, vjust = 0.5, 
                                                  size = 10, color = "black"), 
                      axis.text.y = element_text (angle = 0, hjust = 1, vjust = 0.5, 
                                                  size = 10, color = "black"),
                      axis.title.y = element_text(size = 12),
                      axis.ticks = element_line(size = 0.2),
                      axis.ticks.length = unit(0.7, 'mm'),
                      axis.line = element_line(size = 0.3),
                      legend.text = element_text(size = 9),
                      legend.title = element_text(size = 11),
  ) +
  scale_x_discrete(limits = c("circumacicular", "acicula", "chaeta",  
                              "muscle", "circumchaetal", "hemichaetal",
                              "epithelia_cell", "ciliated cell"))


w_degree_plot <- as_tibble(desmo_conn_graph %>%  
            as_tbl_graph()) %>%
  filter(class != "other") %>%
  filter(class != "basal lamina") %>%
  select(class, weighted_degree) %>%
  ggplot(aes(x = class, y = weighted_degree)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.25, size = 0.6, aes(color = weighted_degree, alpha = weighted_degree/84)) +
  scale_colour_gradient(low = Reds[1], high = Reds[9], 
                        guide = "colorbar", 
                        guide_legend('weighted degree')) +
  labs(y = 'weighted degree', x = "") +
  theme_half_open() + 
  guides(alpha = "none", color = "colorbar") +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text (angle = 90, hjust = 1, vjust = 0.5, 
                                    size = 10, color = "black"), 
        axis.text.y = element_text (angle = 0, hjust = 1, vjust = 0.5, 
                                    size = 10, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.7, 'mm'),
        axis.line = element_line(size = 0.3),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(4, "mm"),
        legend.position = 'top'
  ) +
  scale_x_discrete(limits = c("circumacicular", "acicula", "chaeta",  
                              "muscle", "circumchaetal", "hemichaetal",
                              "epithelia_cell", "ciliated cell"),
                   labels = c("circumacicular", "aciculae", "chaetae",  
                              "muscles", "circumchaetal", "hemichaetal",
                              "epidermal cells", "ciliated cells"))

ggsave("pictures/Figure6_w_degree_plot.png", w_degree_plot,
       width = 1200, units = "px",  height = 1200, limitsize = TRUE)

}

# read and analyse grouped graph ------------------------------------------

#read the saved visNetwork file from source_data/
desmo_grouped_graph <- readRDS("source_data/Figure3_figure_supplement1_source_data3.rds")

# generate and plot subgraphs ---------------------------------------------------------------

{

#define color list
cols <- c(M_Winding_Col, Okabe_Ito, Reds, brewer12, Tol_muted)
group_colors <- sample(cols, length(V(desmo_grouped_graph)), replace = TRUE)
V(desmo_grouped_graph)$color <- group_colors

#get subgraph of a single node defined by class
subgraph <- function(class_id) {desmo_grouped_graph %>% 
  as_tbl_graph() %>%
  convert(to_local_neighborhood,
          node = which(.N()$class == class_id),
          order = 1,
          mode = "all")
}

#derive subgraphs
subgraph_CA <- subgraph('circumacicular') %>%
  toVisNetworkData()
subgraph_A <- subgraph('acicula') %>%
  toVisNetworkData()
subgraph_C <- subgraph('chaeta') %>%
  toVisNetworkData()
subgraph_CC <- subgraph('circumchaetal') %>%
  toVisNetworkData()

#write edge weight to value column (for vis)
subgraph_CA$edges$value <- subgraph_CA$edges$weight
subgraph_CC$edges$value <- subgraph_CC$edges$weight
subgraph_C$edges$value <- subgraph_C$edges$weight
subgraph_A$edges$value <- subgraph_A$edges$weight

#add names to label column
subgraph_CA$nodes$label <- subgraph_CA$nodes$CATMAID_name
subgraph_CC$nodes$label <- subgraph_CC$nodes$CATMAID_name
subgraph_C$nodes$label <- subgraph_C$nodes$CATMAID_name
subgraph_A$nodes$label <- subgraph_A$nodes$CATMAID_name

#add class as group label
subgraph_CA$nodes$group <- subgraph_CA$nodes$class
subgraph_CC$nodes$group <- subgraph_CC$nodes$class
subgraph_C$nodes$group <- subgraph_C$nodes$class
subgraph_A$nodes$group <- subgraph_A$nodes$class

#define function for visNet graph visualisation
visNet <- function(nodes, edges) {visNetwork(nodes, 
                     edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE) %>%
  visPhysics(solver = "forceAtlas2Based", 
             forceAtlas2Based = list(gravitationalConstant = -150000,
                                     centralGravity = 10,
                                     springConstant = 0.001,
                                     avoidOverlap = 1)) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=1, max=25),
           color = list(color = "#454545", opacity = 0.2)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           scaling=list(min=20, max=50),
           font = list(size = 35, 
                       color = "black", 
                       background = "#DEDEDE")) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE)  %>%
  visGroups(groupname = "muscle", color='grey', shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "epithelia_cell", color="#0072B2", shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "circumchaetal", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "circumacicular", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "acicula", color="#D55E00", shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "chaeta", color="#D55E00", shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "hemichaetal", color="#56B4E9", shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "ciliated cell", color="#E69F00", shape = "dot", 
            opacity=1) 
}

visNet(subgraph_CC$nodes, subgraph_CC$edges)
visNet(subgraph_CA$nodes, subgraph_CA$edges)
visNet(subgraph_A$nodes, subgraph_A$edges)
visNet(subgraph_C$nodes, subgraph_C$edges)

#define graph save function
saveVisIgraph <- function(visNetname, graphname) {
#save as html
saveNetwork(visNetname, paste("pictures/Figure6_graph", graphname, ".html", sep = ""), 
            selfcontained = TRUE)
#create png snapshot
webshot2::webshot(url=paste("pictures/Figure6_graph", graphname, ".html", sep = ""),
                  file=paste("pictures/Figure6_graph", graphname, ".png", sep = ""),
                  vwidth = 1500, vheight = 1500, #define the size of the browser window
                  cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}

visCC <- visNet(subgraph_CC$nodes, subgraph_CC$edges)
saveVisIgraph(visCC, 'circumchae')
visCA <- visNet(subgraph_CA$nodes, subgraph_CA$edges)
saveVisIgraph(visCA, 'circumac')
visC <- visNet(subgraph_C$nodes, subgraph_C$edges)
saveVisIgraph(visC, 'acic')
visA <- visNet(subgraph_A$nodes, subgraph_A$edges)
saveVisIgraph(visA, 'chae')

}


# OUTLINE MUS anatomical pics ---------------------------------------------
#read chaetae and aciculae for sg2 only
{
chaeta_sg2_skids <- skids_by_2annotations("chaeta", "segment_2")
chaeta_sg2 <- nlapply(read.neurons.catmaid(chaeta_sg2_skids, pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
acic_sg2_skids <- skids_by_2annotations("acicula", "segment_2")
acicula_sg2 <- nlapply(read.neurons.catmaid(acic_sg2_skids, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
CC_sg2_skids <- skids_by_2annotations("circumchaetal", "segment_2")
CC_sg2 <- nlapply(read.neurons.catmaid(CC_sg2_skids, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
CA_sg2_skids <- skids_by_2annotations("circumacicular", "segment_2")
CA_sg2 <- nlapply(read.neurons.catmaid(CA_sg2_skids, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
}

#parapodial background plotting function
plot_background_sg2_parapod <- function(){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.07,
         col="#E2E2E2")
  plot3d(acicula_sg2, soma=T, lwd=4,
         add=T, alpha=0.5,
         col="black")
  plot3d(chaeta_sg2, soma=T, lwd=2,
         add=T, alpha=1,
         col="grey80")
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.55)
  nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 1))
  #z-axis clip
  clipplanes3d(0, 0, -1, 133000)
  #z-axis clip from top
  clipplanes3d(0, 0, 1, -66000)
  #y-axis clip
  clipplanes3d(1, 0, 0.001, -57000)
  #x-axis clip
  clipplanes3d(0.1, 1, 0, -57000)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  #rotation
  nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 1))
  play3d( spin3d( axis = c(1, 0, 0), rpm = 15), duration = 1)
  play3d( spin3d( axis = c(0, 1, 0), rpm = 1), duration = 1)
  play3d( spin3d( axis = c(0, 0, 1), rpm = 1), duration = 1)
  
}

#read all outlines
OUTLINES <-   nlapply(read.neurons.catmaid("^outline$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))

#get the skids and make a tibble of skids and name
OUTLINE_skids <- names(OUTLINES)
OUTLINE_names <- catmaid_get_neuronnames(OUTLINE_skids, pid = 11)
OUTLINE_skids <- attributes(OUTLINE_names) #skids are the attributes of the names
OUTLINE_names <- unname(OUTLINE_names)

plot_background_sg2_parapod()
circumac_color <- tibble(subgraph_CA$nodes) %>%
  filter(class == 'circumacicular') %>%
  select(color) %>%
  pull()
plot3d(CA_sg2, soma=T, lwd=4, add=T, alpha=1, col=circumac_color)

#read the MUS OUTLINES for the circumacicular subgraph and plot in the color same as the graph
for (i in 1:length(subgraph_CA$nodes$label)) {
  color <- subgraph_CA$nodes$color[i]
  #only names in sg2
  name_to_filter <- paste(subgraph_CA$nodes$CATMAID_name[i], "_sg2", sep = "")
  MUS_skids <- tibble(skids = unlist(OUTLINE_skids), names = OUTLINE_names) %>%
   filter(grepl(name_to_filter, names)) %>%
    select(skids) %>%
    pull()
  if (length(MUS_skids) == 0) {next} #some names may not return skids, skip these
  MUS_OUTLINE <- read.neurons.catmaid(MUS_skids, pid=11)
  plot3d(MUS_OUTLINE, soma=T, lwd=10, add=T, alpha=0.6, col=color)
}

#save pic
rgl.snapshot("pictures/Fig6_sg2para_MUS_CA_1.png")
close3d()

plot_background_sg2_parapod()
circumchaetal_color <- tibble(subgraph_CC$nodes) %>%
  filter(class == 'circumchaetal') %>%
  select(color) %>%
  pull()
plot3d(CC_sg2, soma=T, lwd=4, add=T, alpha=0.5, col=circumchaetal_color)

#read the MUS OUTLINES for the circumachaetal subgraph and plot in the color same as the graph
for (i in 1:length(subgraph_CC$nodes$label)) {
  color <- subgraph_CC$nodes$color[i]
  #only names in sg2
  name_to_filter <- paste(subgraph_CC$nodes$CATMAID_name[i], "_sg2", sep = "")
  MUS_skids <- tibble(skids = unlist(OUTLINE_skids), names = OUTLINE_names) %>%
    filter(grepl(name_to_filter, names)) %>%
    select(skids) %>%
    pull()
  if (length(MUS_skids) == 0) {next} #some names may not return skids, skip these
  MUS_OUTLINE <- read.neurons.catmaid(MUS_skids, pid=11)
  plot3d(MUS_OUTLINE, soma=T, lwd=10, add=T, alpha=0.6, col=color)
}

#save pic
rgl.snapshot("pictures/Fig6_sg2para_MUS_CC_1.png")
close3d()

# assemble figure ---------------------------------------------------------

{
panel_A <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_desmo_connectome_highest_weight_ventral.png")) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
           x = 0.83, y = 0.14, size = 8) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.82),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.82,
                   xend = 0.1,
                   yend = 0.9),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.93, size = 7) +
  draw_label("p", x = 0.1, y = 0.80, size = 7) 

panel_B <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_desmo_connectome_highest_weight_right.png")) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.25,
                   yend = 0.9),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.25,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.9),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("v", x = 0.06, y = 0.9, size = 7) +
  draw_label("d", x = 0.29, y = 0.9, size = 7) 

panel_degree_gr <- ggdraw() + draw_image(readPNG("pictures/Fig6_desmo_connectome_degree.png")) 

degree_scale <- ggdraw() + draw_image(readPNG("pictures/Figure6_degree_scale.png"))

panel_w_degree <- ggdraw() + draw_image(readPNG("pictures/Figure6_w_degree_plot.png"))

panel_CC <- ggdraw() + draw_image(readPNG("pictures/Figure6_graphcircumchae.png")) + 
  draw_label("circumchaetal & partners", x = 0.5, y = 0.97, size = 8)
panel_CA <- ggdraw() + draw_image(readPNG("pictures/Figure6_graphcircumac.png")) + 
  draw_label("circumacicular & partners", x = 0.5, y = 0.97, size = 8)
panel_Ac <- ggdraw() + draw_image(readPNG("pictures/Figure6_graphacic.png")) + 
  draw_label("aciculae & partners", x = 0.5, y = 0.97, size = 8)
panel_Ch <- ggdraw() + draw_image(readPNG("pictures/Figure6_graphchae.png")) + 
  draw_label("chaetae & partners", x = 0.5, y = 0.97, size = 8)
  

panel_Circumacic_partners <- ggdraw() +
  draw_label(paste("25 ", "\u00B5", "m", sep = ""), 
             x = 0.83, y = 0.14, size = 9) + 
  draw_image(readPNG("pictures/Fig6_sg2para_MUS_CA_1.png")) + 
  draw_label("circumacicular & partners", x = 0.5, y = 0.97, size = 8)  + 
  draw_label("neuropodium", x = 0.82, y = 0.88, size = 7)  + 
  draw_label("notopodium", x = 0.8, y = 0.1, size = 7)  +
  geom_segment(aes(x = 0.1,
                   y = 0.8,
                   xend = 0.1,
                   yend = 0.72),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.72,
                   xend = 0.1,
                   yend = 0.8),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.83, size = 7) +
  draw_label("p", x = 0.1, y = 0.70, size = 7) 

panel_Circumchaetal_partners <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_sg2para_MUS_CC_1.png")) + 
  draw_label("circumchaetal & partners", x = 0.5, y = 0.97, size = 8) + 
  draw_label("neuropodium", x = 0.82, y = 0.88, size = 7)  + 
  draw_label("notopodium", x = 0.8, y = 0.1, size = 7) +
  geom_segment(aes(x = 0.1,
                   y = 0.8,
                   xend = 0.1,
                   yend = 0.72),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.72,
                   xend = 0.1,
                   yend = 0.8),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.83, size = 7) +
  draw_label("p", x = 0.1, y = 0.70, size = 7) 

layout <- "
AAAAAAABBBBBBBCCCCCCCCCCDDDDDDDD
EEEEEEEEFFFFFFFFGGGGGGGGHHHHHHHH
"

Figure6 <- panel_A + panel_B + panel_degree_gr + panel_w_degree +
  panel_Circumacic_partners + panel_CA + panel_Circumchaetal_partners  + panel_CC +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure6.pdf", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2400, height = 1400)

}

ggsave("figures/Figure6.png", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2400, height = 1400, bg = 'white')
