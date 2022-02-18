
#load packages
{
library(igraph)
library(visNetwork)
library(ggplot2)
library(ggraph)
library(rgexf)
library(RColorBrewer)
}
 

#read adjacency matrix
desmo_conn <- read.csv("data/adjacency_matrix_desmosomal_connectome_CATMAID.csv", sep=",")
desmo_conn <- as.matrix(desmo_conn)

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  desmo_conn[,2:dim(desmo_conn)[2]],
  mode = c("undirected"), weighted = TRUE, diag = TRUE)

#create a new graph only from the largest connected component
#check connected componenets
cl <- components(Conn_graph)
g_graph_largest = induced_subgraph(g, which(cl$membership == 1))
g_graph_largest

#calculate node weights
weights=strength(g_graph_largest, vids = V(g_graph_largest), mode = c("all"),
                 loops = TRUE)
degree=degree(g_graph_largest, v = V(g_graph_largest), mode = c("all"), loops = TRUE, normalized = FALSE)

# passing some info  
V(g_graph_largest)$value <- degree
E(g_graph_largest)$value <- E(g_graph_largest)$weight



# Serialise the graph to Gephi format
#gexf_data <- rgexf::igraph.to.gexf(g_graph_largest)
# Write the network into a gexf (Gephi) file
#write.gexf(gexf_data, output = "data/g_graph_largest.gexf")
#do force field layout and Leiden modularity in Gephi

#reimport gephi file after layout in gephi
#gexf_data <- read.gexf("data/g_graph_largest_layout.gexf")

#convert to igraph
#g_graph_largest <- rgexf::gexf.to.igraph(gexf_data)

# Compute Leiden communitites from the graph data
com <- cluster_leiden(g_graph_largest,
                      objective_function = "modularity",
                      weights = NULL,
                      resolution_parameter = 0.3,
                      beta = 0.01,
                      n_iterations = 100)
#number of communities
N_communities <- max(com$membership)
N_communities


#generate coluors for the different communities
col=list()
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
col[1:7] <- Okabe_Ito[1:7]
col[8:max(com$membership)]=rev(brewer.pal(max(com$membership)-7, "Paired"))

#get coordinates from imported gephi file
coords <- as.matrix(data.frame(x = gexf_data$nodesVizAtt$position$x, y = gexf_data$nodesVizAtt$position$y, stringsAsFactors = FALSE))

## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(g_graph_largest)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight

#assign node membership
Conn_graph.visn$nodes$group <- com$membership
#colour by group
Conn_graph.visn$nodes$color <- unlist(lapply(Conn_graph.visn$nodes$group, function(x) col[Conn_graph.visn$nodes$group[x]]))

Conn_graph.visn$nodes$size <- sqrt(V(g_graph_largest)$value)*10


#plot the same graph by using the coordinates from the gephi layout
visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 42) %>%
  visNodes(borderWidth=0.3, 
           color = list(background=Conn_graph.visn$nodes$color, border='black'),
           opacity=1,
           shape='dot',
           physics=TRUE,
           font=list(color='black', size=0),
           scaling = list(label=list(enabled=TRUE, min=16, max=28))) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
           scaling=list(min=2, max=8),
           color = list(inherit=TRUE, opacity=0.4)) %>%
  visPhysics(solver = "forceAtlas2Based",
           forceAtlas2Based = list(gravitationalConstant = -1000, 
                     theta=0.1,
                     springLength=30,
                     avoidOverlap=1,
                     centralGravity=0.5,
                     springConstant=0.5,
                     damping=1))

visNet
