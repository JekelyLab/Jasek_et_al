require(igraph)
library(natverse)
library(visNetwork)
library(ggplot2)
library(ggraph)

 
conn <- source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))

#read adjacency matrix

desmo_conn <- read.csv("/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Muscles/code/Jasek_et_al/data/adjacency_matrix_desmosomal_connectome_CATMAID.csv",
         sep=",")
desmo_conn <- as.matrix(desmo_conn)


#alernatively, retrieve neurons from catmaid
connectome_IN = read.neurons.catmaid("^eye_circuit$", pid=11, fetch.annotations = T, conn=conn_http1)

skids <- attr(connectome_IN,"df")$skid
#or
skids <- names(connectome_IN)
N_neurons <- length(skids)
connectivity_IN <- catmaid_get_connectors_between(pre=skids, 
                                                  post=skids, pid=11)
connectivity_IN
N_neurons
# use table() to cross-tabulate number of connections between skids
connectivity_IN_table = table(connectivity_IN[,1:2])
connectivity_IN_table[1:10,1:10]
dim(connectivity_IN_table)
# create empty connectivity matrix
conn_mat = matrix(0,nrow=N_neurons,ncol=N_neurons)
# name rows and cols for skids
colnames(conn_mat) = skids
rownames(conn_mat) = skids

# populate connectivity matrix based on np_conn_table
for (x in rownames(connectivity_IN_table)){
  for (y in colnames(connectivity_IN_table)){
    conn_mat[x,y] = connectivity_IN_table[x,y]
  }
}

#change skids to neuron names
neuron_names <- catmaid_get_neuronnames(skids, pid = 11, conn = NULL)
colnames(conn_mat) = neuron_names
rownames(conn_mat) = neuron_names

conn_mat
#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  conn_mat, mode = c("directed"),
  weighted = TRUE, diag = TRUE)
#or
Conn_graph <- graph_from_adjacency_matrix(
  desmo_conn[,2:dim(desmo_conn)[2]],
  mode = c("undirected"), weighted = TRUE, diag = TRUE)

#calculate node weights
weights=strength(Conn_graph, vids = V(Conn_graph), mode = c("all"),
                 loops = TRUE)
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

# passing some info  
g <- Conn_graph

#V(g)$color <- c("green")
V(g)$value <- degree*5
V(g)$label.cex = 1
V(g)$label.color = "black"
E(g)$value <- E(g)$weight/10
E(g)[ E(g)$weight > 5 ]$color <- "grey50"

#create a new graph only from the largest connected component
#check connected componenets
cl <- components(g)
g_graph_largest = induced_subgraph(g, which(cl$membership == 1))
g_graph_largest

library(rgexf)
# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(g_graph_largest)
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "g_graph_largest.gexf")
#do force field layout and Leiden modularity in Gephi

#reimport gephi file
gexf_data <- read.gexf("Full_desmosomal_graph.gexf")
#plot as interactive Gephi javascript page
plot(gexf_data)
#convert to igraph
g_graph_largest <- rgexf::gexf.to.igraph(gexf_data)

visIgraph(g_graph_largest)

# Compute clusters of similarity from the graph data
com <- cluster_louvain(g_graph_largest)
#number of communities
N_communities <- max(com$membership)
N_communities
col=hcl.colors(max(com$membership), "Set 2")

# If vertices are in the same cluster, give them the same colour
V(g_graph_largest)$color.border <- "black"
V(g_graph_largest)$color.background <- col[com$membership]
#define font size individually for nodes
V(g_graph_largest)$label.cex <- V(g_graph_largest)$value/30
E(g_graph_largest)$weight
#change vertex shape based on type
#V(g_graph_largest)[V(g_graph_largest)$type == 1]$shape <- "square"
#V(g_graph_largest)[V(g_graph_largest)$type == 0]$shape <- "circle"

#color edges based on vertex color
edge.start <- ends(g_graph_largest, es=E(g_graph_largest),names=F)[,1]
ec <- V(g_graph_largest)$color.background[edge.start]
E(g_graph_largest)$color <- ec


gexf_data$nodesVizAtt$position$x
 
# The DrL layout algorthm uses a random number generator. Setting
# the seed to the same number produces the same layout each time.
{set.seed(42)
  # Comput the DrL graph layout
  l = layout_with_drl(g_graph_largest, options=list(init.attraction=0, 
                                                    simmer.attraction=0,expansion.attraction=0,
                                                    liquid.attraction=0,cooldown.attraction=0,
                                                    cooldown.temperature=0,
                                                    simmer.temperature=600, simmer.iterations=10000
                                                    )) 
}

#define the coordinates based on the drl layout
coords <- as.matrix(data.frame(x = l[,1], y = l[,2], stringsAsFactors = FALSE))
#alternatively, get coordinates from imported gephi file
coords <- as.matrix(data.frame(x = gexf_data$nodesVizAtt$position$x, y = gexf_data$nodesVizAtt$position$y, stringsAsFactors = FALSE))

ec
#V(g_graph_largest)$value <- sqrt(V(g_graph_largest)$value)
# Plot the graph on-screen
plot.igraph(g_graph_largest, layout=coords, frame=FALSE, vertex.size=weights/20, vertex.label=NA, 
     vertex.frame.color='grey25', vertex.label.family = 'sans', edge.width=E(g_graph_largest)$weight/5, 
     edge.arrow.size=E(g_graph_largest)$weight/4, edge.curved=0.1, 
     add=FALSE, edge.color=ec
     #mark.groups = c(5:12),
     #mark.border = 'red',
     #mark.col = 'red',
)




#define node size, shape and color
V(g_graph_largest)$shape <- "circle"
V(g_graph_largest)$value <- degree*50
#define font size individually for nodes
V(g_graph_largest)$label.cex <- V(g_graph_largest)$value/2240
#visualise with visNetwork
visIgraph(g_graph_largest) %>% 
#   visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 42) %>%
#   visPhysics(solver = "forceAtlas2Based",
#             forceAtlas2Based = list(gravitationalConstant = -100, avoidOverlap=1)) %>%
   visEdges(smooth = list(type = 'curvedCW', roundness=0), color='grey') %>%
   visNodes(borderWidth=0.3, color = list(highlight = "yellow")) %>%
   visOptions(highlightNearest = TRUE) %>%
   visInteraction(navigationButtons = TRUE,
                 dragNodes = TRUE, dragView = TRUE,
                 zoomView = TRUE) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords)
  
  



