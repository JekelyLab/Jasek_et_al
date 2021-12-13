require(igraph)
library(visNetwork)
library(ggplot2)
library(ggraph)
library(rgexf)

#reimport gephi file
gexf_data <- read.gexf("data/Full_desmosomal_graph.gexf")

# Serialise igraph formatted graph to Gephi format
#gexf_data <- rgexf::igraph.to.gexf(g_igraph)
# Write the network into a gexf (Gephi) file
#write.gexf(gexf_data, output = "g_igraph.gexf")
#do force field layout and Leiden modularity in Gephi

#plot as interactive Gephi javascript page
#can be opened in browser by clicking "show in new window" in the viewer tab
plot(gexf_data)

#convert to igraph
g_igraph <- rgexf::gexf.to.igraph(gexf_data)

#plotting without customization
visIgraph(g_igraph)

#get coordinates from imported gephi file
coords <- as.matrix(data.frame(x = gexf_data$nodesVizAtt$position$x, y = gexf_data$nodesVizAtt$position$y, stringsAsFactors = FALSE))

#calculate node weighted degree
degree=degree(g_igraph, v = V(g_igraph), mode = c("all"), loops = TRUE, normalized = FALSE)


#define node size, shape and color
V(g_igraph)$shape <- "dot"
V(g_igraph)$value <- degree*15
#define font size individually for nodes
V(g_igraph)$label.cex <- V(g_igraph)$value/200

#add new column with 'value' corresponding to edge weight, 'value' is read by visIgraph for edge weight
E(g_igraph)$value <- E(g_igraph)$weight

#define background and border color separately
V(g_igraph)$color.background <- V(g_igraph)$color
V(g_igraph)$color.border <- "black"

#color edges based on vertex color.background
edge.start <- ends(g_igraph, es=E(g_igraph),names=F)[,1]
ec <- V(g_igraph)$color[edge.start]
E(g_igraph)$color <- ec




#customized visualisation with visNetwork
#can be opened in browser by clicking "show in new window" in the viewer tab
visIgraph(g_igraph) %>% 
   visEdges(smooth = list(type = 'curvedCW', roundness=0), color=ec) %>%
   visNodes(borderWidth=0.2, labelHighlightBold = TRUE) %>%
   visOptions(highlightNearest = TRUE) %>%
   visInteraction(navigationButtons = TRUE,
                 dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords)
  
  



