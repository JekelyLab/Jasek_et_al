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

# plot graph with nodes colored by weighted degree ------------------------

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

visNet_degree <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=5, max=30),
           color = list(source='both', opacity=0.15)) %>%
  visNodes(borderWidth=0.3, 
           color = list(border='black'),
           opacity = 1, 
           scaling=list(min=1, max=55),
           font = list(size = 0)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',
                                     labelOnly=FALSE), 
             width = 1500, height = 1500, autoResize = FALSE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE, hover=TRUE,
                 multiselect=TRUE)

visNet_degree

#make scale of degree
png("pictures/Figure6_degree_scale.png",
    units = c("px"), width = 440, height = 800, bg = 'white')
dotchart(rep(1,9), bg =Reds[1:9], c('1','','','','','','','','84'), 
         frame.plot = FALSE, xaxt = 'n',
         xlim = c(0,4), offset = 2, pt.cex = 12, cex = 4,
         lcolor = 'white', ylab = 'weighted degree')
dev.off()

#save as html
saveNetwork(visNet_degree, "pictures/Fig6_desmo_connectome_degree.html", selfcontained = TRUE)
#create png snapshot
webshot::webshot(url="pictures/Fig6_desmo_connectome_degree.html",
                 file="pictures/Fig6_desmo_connectome_degree.png",
                 vwidth = 1500, vheight = 1500, #define the size of the browser window
                 cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)


# make plot of node weighted degree and degree  ----------------------------

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
  geom_jitter(width = 0.2, size = 0.6, aes(color = weighted_degree, alpha = weighted_degree/84)) +
  scale_colour_gradient(low = Reds[1], high = Reds[9], 
                        guide = "colorbar", 
                        guide_legend('weighted degree')) +
  labs(y = 'weighted degree', x = 'cell class') +
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

# assemble figure ---------------------------------------------------------

panel_A <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_desmo_connectome_highest_weight_ventral.png")) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
           x = 0.83, y = 0.14, size = 9) +
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

panel_B <- ggdraw() + 
  draw_image(readPNG("pictures/Fig6_desmo_connectome_highest_weight_right.png")) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.22,
                   yend = 0.9),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.22,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.9),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("v", x = 0.07, y = 0.9, size = 8) +
  draw_label("d", x = 0.25, y = 0.9, size = 8) 

img_degree_gr <- image_read("pictures/Fig6_desmo_connectome_degree.png")  %>%
  image_rotate(90)
panel_degree_gr <- ggdraw() + draw_image(img_degree_gr) 

degree_scale <- ggdraw() + draw_image(readPNG("pictures/Figure6_degree_scale.png"))
# Add scale as inset element
panel_C <- panel_degree_gr + 
  inset_element(degree_scale, 0.72, 0.72, 1, 1, ignore_tag = TRUE)


layout <- "
ABCD
"

Figure6 <- panel_A + panel_B + panel_C + w_degree_plot +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure6.pdf", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2400, height = 800)


ggsave("figures/Figure6.png", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2600, height = 1500, bg = 'white')
