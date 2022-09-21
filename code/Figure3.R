#Code to generate Figure 3 showing the different Leiden modules of the desmosomal connectome graph from Jasek et al 2021
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

#read the saved igraph format graph file from supplements/
#this file was generated with the desmo_connectome_graph.R script
desmo_conn_graph <- readRDS("source_data/Figure3_source_data1.rds")

#read the saved visNetwork file from supplements/
conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")

# plot graph with coordinates from gephi ----------------------------------

{
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)

conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$weight
conn_graph.visn$nodes$label <- conn_graph.visn$nodes$CATMAID_name

#save for supplement
visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=2, max=15),
             color = list(inherit=TRUE, opacity=0.5)) %>%
    visNodes(borderWidth=0.3, 
             color = list(border='black'),
             opacity = 1, 
             scaling=list(min=15, max=50),
             font = list(size = 15)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',
                                       labelOnly=FALSE), 
               width = 1500, 
               height = 1500,
               autoResize = TRUE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE)

visNet
#save as html
saveNetwork(visNet, "source_data/Figure3_source_data3.html", selfcontained = TRUE)

#save for figure without node labels
visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=2, max=15),
           color = list(inherit=TRUE, opacity=0.5)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           opacity = 1, 
           scaling=list(min=15, max=50),
           font = list(size = 0)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, 
             height = 1500,
             autoResize = TRUE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE)

#save as html
saveNetwork(visNet, "pictures/Full_desmo_connectome_modules_webshot.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Full_desmo_connectome_modules_webshot.html",
                  file="pictures/Full_desmo_connectome_modules_webshot.png",
                  vwidth = 1500,
                  vheight = 1500, #define the size of the browser window
                  cliprect = c(50, 60, 1500, 1500), 
                  zoom=5, 
                  delay = 2
                  )
}

# plot neurons by colours matching network pic -----------------------------
{
#get names and colours and make df
names <- sapply(conn_graph.visn$nodes$id, function(x) x)
colors <- sapply(conn_graph.visn$nodes$color, function(x) x)
module <- sapply(conn_graph.visn$nodes$partition, function(x) x)
df <- data.frame(names, colors, module)
  
#get skid for module 4, 2nd cell
df[df$module==4,][2,][1]
  
#iterate through the neurons and plot them in the colour that was used in the network plot
{
for (i in 1:13){
      print(i)
      plot_background_ventral()
      for (j in 1:nrow(df[df$module==i,])){
        #read the skeleton based on the skid
        skeleton <- nlapply(read.neurons.catmaid(df[df$module==i,][j,][[1]], pid=11),
                            function(x) smooth_neuron(x, sigma=6000))
        #plot skeleton
        plot3d(skeleton, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
               rev = FALSE, fixup = F, add=T, alpha=1, col=df[df$module==i,][j,2])
        }
      filename = paste("pictures/desmo_conn_module_", i, ".png", sep = "")
      rgl.snapshot(filename)
      close3d()
     }
    }
}


# graph statistics --------------------------------------------------
#assign degree to nodes
V(desmo_conn_graph)$degree <-
degree(
  desmo_conn_graph,
  v = V(desmo_conn_graph),
  mode = "all",
  loops = TRUE,
  normalized = FALSE
)

#assign weighted degree to nodes
V(desmo_conn_graph)$weighted_degree <-
  strength(
    desmo_conn_graph,
    v = V(desmo_conn_graph),
    mode = "all",
    loops = TRUE,
)

desmo_conn_graph.tb <- desmo_conn_graph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate("authority" = centrality_authority(weights = E(desmo_conn_graph)$weight)) %>%
  mutate("pagerank" = centrality_pagerank(weights = E(desmo_conn_graph)$weight, 
                                          directed = TRUE)) %>%
  mutate("closeness" = centrality_closeness(weights = E(desmo_conn_graph)$weight,
                                            mode = "all")) %>%
  mutate("eigen" = centrality_eigen(directed = TRUE)) %>%
  mutate("hub" = centrality_hub(weights = E(desmo_conn_graph)$weight) ) %>%
  mutate("node_eccentricity" = node_eccentricity(mode = "all") ) %>%
  mutate("node_is_cut" = node_is_cut() )  %>%
  mutate("local_triangles" = local_triangles())

#plot node degree distribution
{
plot_degree <- as.data.frame(V(desmo_conn_graph)$degree) %>%
    ggplot(aes(x=V(desmo_conn_graph)$degree)) +
    geom_histogram(binwidth = 0.05) +
    scale_x_log10(breaks = c(1,2,5, 10,20))  +
    labs(x="degree",y="# of nodes",title=" ") +
    theme_minimal_hgrid() +
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=8),
          panel.grid.minor = element_blank() ) 
plot_degree  
}

#plot node degree against eccentricity
{
as_tibble(desmo_conn_graph.tb %>%
    filter(group %in% c("muscle", "acicula", "circumacicular",
                        "chaeta", "circumchaetal", "hemichaetal"))) %>%
    ggplot(aes(x = degree, 
               y = unlist(pagerank), 
               color = group,
               alpha = group,
               shape = group,
               size = node_eccentricity)) +
    geom_point() +
    scale_x_log10()  +
    labs(x="degree",y="pagerank",
         title=" ", color = 'cell class', size = 'eccentricity',
         alpha = 'cell class', shape = 'cell class') +
    theme_minimal_hgrid() +
    scale_color_manual(values = list(muscle = "grey50",
                                     acicula = "#CC79A7", 
                                     circumacicular = "#0072B2", 
                                     chaeta = "#E69F00", 
                                     circumchaetal = 'black',
                                     hemichaetal = 'red'),
                       labels = c("muscle", "acicula", 
                                  "circumacicular", 
                                  "chaeta", 
                                  "circumchaetal", 
                                  "hemichaetal") )  +
    scale_shape_manual(values = list(muscle = 1, 
                                     acicula = 16, 
                                     circumacicular = 15, 
                                     chaeta = 17, 
                                     circumchaetal = 18,
                                     hemichaetal = 19),
                       labels = c("muscle", 
                                  "acicula", 
                                  "circumacicular", 
                                  "chaeta", 
                                  "circumchaetal", 
                                  "hemichaetal")) +
    scale_alpha_manual(values = c(muscle = 0.2, 
                                     acicula = 0.7, 
                                     circumacicular = 0.6, 
                                     chaeta = 0.3, 
                                     circumchaetal = 0.3,
                                     hemichaetal = 0.5),
                       labels = c("muscle", 
                                  "acicula", 
                                  "circumacicular", 
                                  "chaeta", 
                                  "circumchaetal", 
                                  "hemichaetal") ) +
    theme(legend.position='right', 
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text = element_text(size=12),
          legend.text = element_text(size=8),
          legend.title = element_text(size=9),
          legend.key.size = unit(3, "mm"),
          panel.grid.minor = element_blank() ) +
    guides(alpha = 'none')

ggsave("pictures/connectome_weighted_degree_vs_pagerank.png", 
         width=1600, height = 1200, units = "px", bg = "white")
}

#plot distribution of cell classes across partitions 
{

#overwrite group with type_of_cell (class)  
plot_modules <- as_tibble(desmo_conn_graph.tb) %>%
    group_by(partition) %>%         # Specify group indicator
    filter(group %in%  c("muscle", "acicula", "circumacicular",
                         "chaeta", "circumchaetal", "hemichaetal",
                         "ciliated cell",
                         "epithelia_cell") ) %>%
    ggplot(aes(x = partition, fill = group)) +
    geom_bar() +
    scale_fill_manual(values = list(muscle = "grey50",
                                    acicula = "#CC79A7", 
                                    circumacicular = "#0072B2", 
                                    chaeta = Okabe_Ito[1], 
                                    circumchaetal = 'black',
                                    hemichaetal = Okabe_Ito[2],
                                    `ciliated cell` = Okabe_Ito[6],
                      `epithelia_cell` = Okabe_Ito[3]),
                      labels = c("muscle", "acicula", 
                                 "circumacicular", 
                                 "chaeta", 
                                 "circumchaetal",
                                 "hemichaetal",
                                 "ciliated cell", "EC")  ) +
    scale_x_continuous(limits = c(0,13), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
    labs(x = 'module', y = '# of cells', fill = 'cell class') +
    theme_minimal_hgrid()+
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=8),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10),
          legend.key.size = unit(3, "mm"),
          panel.grid.minor = element_blank(),
          legend.position = 'right')

plot_modules

}


# make table about cell statistics -------------------------------------

{
  
#get counts of cell types and nodes/edges
N_mus <- dim(as_tibble(desmo_conn_graph.tb) %>%
    group_by(partition) %>%         # Specify group indicator
    filter(group %in%  c("muscle") ))[1]

N_bl <- dim(as_tibble(desmo_conn_graph.tb) %>%
               group_by(partition) %>%         # Specify group indicator
               filter(group %in%  c("basal lamina") ))[1]
N_bl

N_EC <- dim(as_tibble(desmo_conn_graph.tb) %>%
               group_by(partition) %>%         # Specify group indicator
               filter(group %in%  c("epithelia_cell") ))[1]

N_bl <- dim(as_tibble(desmo_conn_graph.tb) %>%
              group_by(partition) %>%         # Specify group indicator
              filter(group %in%  c("basal lamina") ))[1]

N_ch_ac <- dim(as_tibble(desmo_conn_graph.tb) %>%
      group_by(partition) %>%         # Specify group indicator
      filter(group %in%  c("chaeta", "acicula") ))[1]

N_cil <- dim(as_tibble(desmo_conn_graph.tb) %>%
                 group_by(partition) %>%         # Specify group indicator
                 filter(group %in%  c("ciliated cell") ))[1]

N_circ <- dim(as_tibble(desmo_conn_graph.tb) %>%
      group_by(partition) %>%         # Specify group indicator
      filter(group %in%  c("circumacicular", "circumchaetal") ))[1]

N_hemi <- dim(as_tibble(desmo_conn_graph.tb) %>%
                group_by(partition) %>%         # Specify group indicator
                filter(group %in%  c("hemichaetal") ))[1]

N_node <- length(desmo_conn_graph.tb %>%
      select(name)) 

N_frag <- dim(as_tibble(desmo_conn_graph.tb) %>%
  filter(segment == "fragmentum") )[1]

N_edge <- length(E(as.igraph(desmo_conn_graph.tb)))

N_desmo <- desmo_conn_graph.tb %>%
  activate(edges) %>%
  select(weight) %>%
  pull() %>%
  sum()

table <- plot_ly(
    type = 'table',
    columnwidth = c(10, 7),
    columnorder = c(0, 1),
    header = list(
      values = c("Type","Number"),
      align = c("center", "center"),
      line = list(width = 1, color = 'black'),
      fill = list(color = c("#E69F00", "#0072B2")),
      font = list(family = "Arial", size = 14, color = "white")
    ),
    cells = list(
      values = rbind(c("total nodes in the connectome", 
                       "nodes with soma",
                       "edges", 
                       "in-graph desmosomes",
                       "muscle cells",
                       "epidermal cells", 
                       "aciculo-/ chaetoblasts", 
                       "circumacicular/-chaetal cells",
                       "ciliary band cells",
                       "hemichaetal"), 
                     c(N_node,
                       N_node-N_frag,
                       N_edge,
                       N_desmo,
                       N_mus, 
                       N_EC, 
                       N_ch_ac,
                       N_circ, 
                       N_cil,
                       N_hemi)),
      align = c("center", "center"),
      line = list(color = "black", width = 0.3),
      font = list(family = "Arial", size = 13, color = c("black"))
    ))
  
table

saveNetwork(table, "pictures/desmo_connectome_stats_table.html")
webshot::webshot(url="pictures/desmo_connectome_stats_table.html",
                    file="pictures/desmo_connectome_stats_table.png",
                    vwidth=400, vheight=300, #define the size of the browser window
                    cliprect = c(20,50,350, 220), zoom=10)
  
}

# assemble figure ---------------------------------------------------------

{

#read with image_read (magick package) and rotate
img_conn <- image_rotate(image_read("pictures/Full_desmo_connectome_modules_webshot.png"), 270) %>%
  image_flip()

panel_mod1 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_1.png")) +
    draw_label("1) oblique (sg2)", x = 0.5, y = 0.05, size = 9)
panel_mod2 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_2.png")) +
    draw_label("2) oblique (sg3)", x = 0.5, y = 0.05, size = 9)
panel_mod3 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_3.png")) +
    draw_label("3) ventrolateral (l)", x = 0.5, y = 0.05, size = 9)
panel_mod4 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_4.png")) +
    draw_label("4) parapodial (sg2r)", x = 0.5, y = 0.05, size = 9)
panel_mod5<- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_5.png")) +
    draw_label("5) parapodial (sg1r)", x = 0.5, y = 0.05, size = 9)
panel_mod6 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_6.png")) +
    draw_label("6) sg0 and head", x = 0.5, y = 0.05, size = 9)
panel_mod7 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_7.png")) +
    draw_label("7) dorsolateral (l)", x = 0.5, y = 0.05, size = 9)
panel_mod8 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_8.png")) +
    draw_label("8) parapodial (sg1l)", x = 0.5, y = 0.05, size = 9)
panel_mod9 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_9.png")) +
    draw_label("9) parapodial (sg3r)", x = 0.5, y = 0.05, size = 9)
panel_mod10 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_10.png")) +
  draw_label("10) parapodial (sg3l)", x = 0.5, y = 0.05, size = 9)
panel_mod11 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_11.png")) +
    draw_label("11) dorsolateral (r)", x = 0.5, y = 0.05, size = 9)
panel_mod12 <- ggdraw() + draw_image(readPNG("pictures/desmo_conn_module_12.png")) +
    draw_label("12) ventrolateral (r)", x = 0.5, y = 0.05, size = 9)

panel_conn <- ggdraw() + draw_image(img_conn) +
    draw_label("sg0 and head", x=0.42, y = 0.9, size = 9) +
    draw_label("parapodial (sg1l)", x=0.65, y = 0.72, size = 9) +
    draw_label("parapodial (sg1r)", x=0.27, y = 0.65, size = 9) +
    draw_label("oblique (sg2)", x=0.49, y = 0.58, size = 9) +
    draw_label("oblique (sg3)", x=0.54, y = 0.4, size = 9) +
    draw_label("parapodial (sg2r)", x=0.3, y = 0.48, size = 9) +
    draw_label("ventrolateral (l)", x=0.88, y = 0.78, size = 9) +
    draw_label("ventrolateral (r)", x=0.25, y = 0.8, size = 9) +
    draw_label("dorsolateral (l)", x=0.88, y = 0.42, size = 9)  +
    draw_label("dorsolateral (r)", x=0.17, y = 0.32, size = 9) +
    draw_label("parapodial (sg3l)", x=0.8, y = 0.25, size = 9) +
    draw_label("parapodial (sg3r)", x=0.48, y = 0.22, size = 9) 
  
panel_table <- ggdraw() + draw_image(readPNG("pictures/desmo_connectome_stats_table.png"))
  

#define layout with textual representation
layout_A <- "
#BCD#
#BCD#
#BCD#
HIIIJ
HIIIJ
HIIIJ
KIIIL
KIIIL
KIIIL
MIIIN
MIIIN
MIIIN
#PQR#
#PQR#
#PQR#
"
panel_A <-  panel_mod5 + panel_mod6 + panel_mod8 + 
    panel_mod12 + panel_conn + panel_mod3 +
    panel_mod4 + panel_mod1 + 
    panel_mod11 + panel_mod7 + 
    panel_mod9 + panel_mod2 + panel_mod10 + 
    plot_layout(design = layout_A)
#convert into single panel
panel_A <- ggdraw(panel_A)

layout <- "
AAAAAAAAAAAAAAAAAA
BBBBBBBBBCCCCCDDDD
"
  
Fig_conn <- panel_A + 
  panel_table + plot_modules + plot_degree  +
  plot_layout(design = layout, heights = c(3.8,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("figures/Figure3.pdf", limitsize = FALSE, 
         units = c("px"), Fig_conn, width = 2600, height = 3000)
  
}

ggsave("figures/Figure3.png", limitsize = FALSE, 
       units = c("px"), Fig_conn, width = 2600, height = 3000, bg='white')




