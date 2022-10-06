#code for Figure6 of the Jasek et al 2021 desmosomal connectome paper
#quantification and visualisation of node degrees in the desmosomal connectome 
#Gaspar Jekely

source("code/Packages_and_Connection.R")

#read the saved igraph format graph file from supplements/
desmo_conn_graph <- readRDS("source_data/Figure3_source_data1.rds")

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
plot_background_ventral_no_acic()
par3d(zoom=0.55)
  
#plot skeletons with same alpha for comparison
skeletons <- nlapply(read.neurons.catmaid(node_name_and_degree_high_only$name, pid=11), 
                        function(x) smooth_neuron(x, sigma=6000))
plot3d(skeletons, soma = TRUE, lwd = 1, add=T, alpha=0.6, col=Reds[8])
plot3d(outline, add = TRUE, alpha = 0.017)
plot3d(scalebar_50um_ventral, color = "black", lw = 2)
  
#make snapshot
rgl.snapshot("pictures/Fig6_desmo_connectome_same_color.png")
close3d()

#plot skeletons with transparency inversely proportional to weighted degree
plot_background_ventral_no_acic()
par3d(zoom=0.55)

for (i in 1:length(node_name_and_degree_high_only$name)) {
  skeleton <- nlapply(read.neurons.catmaid(node_name_and_degree_high_only$name[i], pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
  plot3d(skeleton, soma = TRUE, lwd = 1, add=T, 
       # alpha=node_name_and_degree_high_only$weighted_degree[i]/max(node_name_and_degree_high_only$weighted_degree),
       alpha = 0.7, 
       #pick discrete colors based on normalised degree
        col=Reds[round(node_name_and_degree_high_only$weighted_degree[i]/max(node_name_and_degree_high_only$weighted_degree)*10)] 
)
}

plot3d(outline, add = TRUE, alpha = 0.017)

#make snapshot
rgl.snapshot("pictures/Fig6_desmo_connectome_highest_weight_ventral.png")

#lateral view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.61)
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
saveNetwork(visNet_degree, "supplements/Fig6_desmo_connectome_degree.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="supplements/Fig6_desmo_connectome_degree.html",
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
  scale_x_discrete(limits = c("acFC", "acicula", "chaeta",  
                              "muscle", "chaeFC-all",
                              "epithelia_cell", "ciliated cell"),
                   labels = c("acFC", "aciculae", "chaetae",  
                              "muscles", "chaeFC",
                              "epidermal cells", "ciliated cells"))

w_degree_plot <- as_tibble(desmo_conn_graph %>%  
            as_tbl_graph()) %>%
  filter(class != "other") %>%
  filter(class != "basal lamina") %>%
  select(class, weighted_degree) %>%
  ggplot(aes(x = class, y = weighted_degree, alpha = 0)) +
  geom_jitter(width = 0.25, size = 0.6, aes(color = weighted_degree, alpha = weighted_degree/84)) +
  geom_boxplot(outlier.shape = NA) + 
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
  scale_x_discrete(limits = c("acFC", "acicula", "chaeta",  
                              "muscle", "chaeFC-all",
                              "epithelia_cell", "ciliated cell"),
                   labels = c("acFC", "aciculae", "chaetae",  
                              "muscles", "chaeFC",
                              "epidermal cells", "ciliated cells"))

ggsave("pictures/Figure6_w_degree_plot.png", w_degree_plot,
       width = 1200, units = "px",  height = 1200, limitsize = TRUE)

#save source data
w_degree_plot_source_data <- as_tibble(desmo_conn_graph %>%  
                                         as_tbl_graph()) %>%
  filter(class != "other") %>%
  filter(class != "basal lamina") %>%
  select(class, weighted_degree)

write_csv(w_degree_plot_source_data, "source_data/Figure6_source_data_1.txt")
#read_csv("source_data/Figure6_source_data_1.txt")

}

# read and analyse grouped graph ------------------------------------------

#read the saved visNetwork file from source_data/
desmo_grouped_graph <- readRDS("source_data/Figure3_figure_supplement1_source_data3.rds")

# generate and plot subgraphs ---------------------------------------------------------------

{
#define colour list
cols <- c(M_Winding_Col, Okabe_Ito, Reds[2:9], brewer12, Tol_muted)
#make sampling reproducible
set.seed(23)
group_colors <- sample(cols, length(V(desmo_grouped_graph)), replace = TRUE)
V(desmo_grouped_graph)$color <- group_colors

#assign color by type
desmo_grouped_graph <- desmo_grouped_graph %>%
  as_tbl_graph() %>%
  mutate(color = ifelse(class == "acicula", "#454545", color)) %>%
  mutate(color = ifelse(class == "chaeta", "#454545", color)) %>%
  mutate(color = ifelse(class == "acFC", "#222222", color)) %>%
  mutate(color = ifelse(class == "chaeFC", "#222222", color))

#get subgraph of a single node defined by class (not used)
subgraph <- function(class_id) {desmo_grouped_graph %>% 
  as_tbl_graph() %>%
  convert(to_local_neighborhood,
          node = which(.N()$class == class_id),
          order = 1,
          mode = "all")
}

subgraph_i <- function(class_id) {
  #make a list of ego graphs for each node selecting only its order 1 neighbourhood
  subgraphs <- desmo_grouped_graph %>% 
    make_ego_graph(order = 1,
                 nodes = which(V(desmo_grouped_graph)$class == class_id),
                 mode = "all")
  #get all the node names from the ego graphs list
  subset_node_names <- unique(unlist(lapply(subgraphs, function(x) print( V(x)$CATMAID_name))))
  #get the corresponding node ids from the original graph
  subset_node_ids <- unlist(lapply(subset_node_names, function(x) which(V(desmo_grouped_graph)$CATMAID_name == x) ) )
  #generate induced subgraph based on the node ids
  desmo_grouped_graph %>% 
    induced.subgraph(vids = subset_node_ids )
}

#derive subgraphs
subgraph_CA <- subgraph_i('acFC') %>%
  toVisNetworkData()
subgraph_CC <- subgraph_i('chaeFC-all') %>%
  toVisNetworkData()

#write edge weight to value column (for vis)
subgraph_CA$edges$value <- subgraph_CA$edges$weight
subgraph_CC$edges$value <- subgraph_CC$edges$weight

#add names to label column
subgraph_CA$nodes$label <- subgraph_CA$nodes$CATMAID_name
subgraph_CC$nodes$label <- subgraph_CC$nodes$CATMAID_name

#add class as group label
subgraph_CA$nodes$group <- subgraph_CA$nodes$class
subgraph_CC$nodes$group <- subgraph_CC$nodes$class
#check label and add level to nodes for hierarchical layout
subgraph_CC$nodes$label
subgraph_CC$nodes$level <- c(1,5,2,
                             2,4,2,
                             4,2,
                             4,3,3,
                             3,3,3,
                             3,3,3,
                             3,3,3,
                             3,3)
subgraph_CA$nodes$label
subgraph_CA$nodes$level <- c(1,7,2,
                             6,2,4,
                             5,4,3,
                             4,3,3,
                             4,3,3,
                             4,4,3,
                             4,5,4,
                             5,5,4,
                             5,5,5,
                             3,3,3,
                             4)
#define function for visNet graph visualisation
visNet <- function(nodes, edges) {visNetwork(nodes, 
                     edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE) %>%
#  visPhysics(solver = "forceAtlas2Based", 
 #            forceAtlas2Based = list(gravitationalConstant = -150000,
  #                                   centralGravity = 6,
   #                                  springConstant = 0.001,
    #                                 avoidOverlap = 1)) %>%
    visHierarchicalLayout(levelSeparation=330, 
                          nodeSpacing=29,
                          direction='LR',
                          sortMethod='hubsize',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=1, max=25),
           color = list(color = "#454545", opacity = 0.2)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           scaling=list(min=20, max=50),
           font = list(size = 42, 
                       color = "black")) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 2500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE)  %>%
  visGroups(groupname = "muscle", shape = "dot", 
            opacity=0.7) %>%
  visGroups(groupname = "epithelia_cell", shape = "dot", 
            opacity=1) %>%
  visGroups(groupname = "chaeFC", size = 30, shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "chaeFC-all", size = 30, shape = "triangle", 
            opacity=1) %>%
  visGroups(groupname = "acicula", size = 25, shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "chaeta", size = 25, shape = "square", 
            opacity=1) %>%
  visGroups(groupname = "ciliated cell", shape = "dot", 
            opacity=1) 
}

#check graphs in Viewer
visNet(subgraph_CA$nodes, subgraph_CA$edges)
visNet(subgraph_CC$nodes, subgraph_CC$edges)


#define graph save function
saveVisIgraph <- function(visNetname, graphname, top, left, width, height) {
#save as html
saveNetwork(visNetname, paste("supplements/Figure6_graph_", graphname, ".html", sep = ""), 
            selfcontained = TRUE)
#create png snapshot
webshot::webshot(url=paste("supplements/Figure6_graph_", graphname, ".html", sep = ""),
                  file=paste("pictures/Figure6_graph_", graphname, ".png", sep = ""),
                  vwidth = 5000, vheight = 5000, #define the size of the browser window
                  cliprect = c(top, left, width, height), zoom=5, delay = 2)

}

visCA <- visNet(subgraph_CA$nodes, subgraph_CA$edges)
saveVisIgraph(visCA, 'acFC', 50,100, 2350, 1500)
visCC <- visNet(subgraph_CC$nodes, subgraph_CC$edges)
saveVisIgraph(visCC, 'chaeFC', 200,500, 1600, 1300)

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
CC_sg2_skids <- skids_by_2annotations("chaeFC-all", "segment_2")
CC_sg2 <- nlapply(read.neurons.catmaid(CC_sg2_skids, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
CA_sg2_skids <- skids_by_2annotations("acFC", "segment_2")
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
acFC_color <- tibble(subgraph_CA$nodes) %>%
  filter(class == 'acFC') %>%
  select(color) %>%
  pull()
plot3d(CA_sg2, soma=T, lwd=4, add=T, alpha=1, col=acFC_color)
#add parapodial scale bar
scale_bar_10um_parapodium_sg2 <- read.neurons.catmaid("^scale_bar_10um_parapodium_sg2$", pid=11)
plot3d(scale_bar_10um_parapodium_sg2, lwd = 4, add = T, col = 'black')

#read the MUS OUTLINES for the acFC subgraph and plot in the color same as the graph
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
chaeFC_color <- tibble(subgraph_CC$nodes) %>%
  filter(class == 'chaeFC') %>%
  select(color) %>%
  pull()
plot3d(CC_sg2, soma=T, lwd=4, add=T, alpha=0.5, col=chaeFC_color[1])

#read the MUS OUTLINES for the chaeFC subgraph and plot in the color same as the graph
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
  draw_image(readPNG("pictures/Fig6_desmo_connectome_same_color.png")) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
           x = 0.85, y = 0.18, size = 8) +
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
  draw_image(readPNG("pictures/Fig6_desmo_connectome_highest_weight_ventral.png"))

panel_C <- ggdraw() + 
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

#degree_scale <- ggdraw() + draw_image(readPNG("pictures/Figure6_degree_scale.png"))

panel_w_degree <- ggdraw() + draw_image(readPNG("pictures/Figure6_w_degree_plot.png"))

panel_CA <- ggdraw() + draw_image(readPNG("pictures/Figure6_graph_acFC.png")) + 
  draw_label("acFC & partners", x = 0.5, y = 0.97, size = 8)  + 
  draw_label("# of desmosomes", x = 0.88, y = 0.82, size = 5) +
  draw_label("1", x = 0.88, y = 0.78, size = 5) + 
  draw_label(max(subgraph_CA$nodes$weighted_degree), x = 0.88, y = 0.74, size = 5) +
  draw_line(x = c(0.9, 0.95), y = c(0.78, 0.78), size = 0.1, color = 'grey') +
  draw_line(x = c(0.9, 0.95), y = c(0.74, 0.74), size = 1.5, color = 'grey')

panel_CC <- ggdraw() + draw_image(readPNG("pictures/Figure6_graph_chaeFC.png")) + 
  draw_label("chaeFC & partners", x = 0.5, y = 0.97, size = 8) + 
  draw_label("# of desmosomes", x = 0.88, y = 0.82, size = 5) +
  draw_label("1", x = 0.84, y = 0.78, size = 5) + 
  draw_label(max(subgraph_CC$nodes$weighted_degree), x = 0.84, y = 0.74, size = 5) +
  draw_line(x = c(0.88, 0.95), y = c(0.78, 0.78), size = 0.1, color = 'grey') +
  draw_line(x = c(0.88, 0.95), y = c(0.74, 0.74), size = 1.1, color = 'grey')

panel_acFC_partners <- ggdraw() +
  draw_image(readPNG("pictures/Fig6_sg2para_MUS_CA_1.png")) + 
  draw_label(paste("10 ", "\u00B5", "m", sep = ""), 
             x = 0.2, y = 0.2, size = 9) + 
  draw_label("acFC & partners", x = 0.5, y = 0.97, size = 8)  + 
  draw_label("notopodium", x = 0.82, y = 0.88, size = 7)  + 
  draw_label("neuropodium", x = 0.8, y = 0.1, size = 7)  +
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
  draw_label("d", x = 0.1, y = 0.83, size = 7) +
  draw_label("v", x = 0.1, y = 0.70, size = 7) 

panel_chaeFC_partners <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_sg2para_MUS_CC_1.png")) + 
  draw_label("chaeFC & partners", x = 0.5, y = 0.97, size = 8) + 
  draw_label("notopodium", x = 0.82, y = 0.88, size = 7)  + 
  draw_label("neuropodium", x = 0.8, y = 0.1, size = 7) +
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
  draw_label("d", x = 0.1, y = 0.83, size = 7) +
  draw_label("v", x = 0.1, y = 0.70, size = 7) 

layout <- "
AAAAAABBBBBBCCCCCCDDDDDDDEEEEEEE
FFFFFFFGGGGGGGGGGGHHHHHHHIIIIIII
"

Figure6 <- panel_A + panel_B + panel_C + panel_degree_gr + panel_w_degree +
  panel_acFC_partners + panel_CA + panel_chaeFC_partners  + panel_CC +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure6.pdf", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2600, height = 1400)

}

ggsave("figures/Figure6.png", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2600, height = 1400, bg = 'white')
