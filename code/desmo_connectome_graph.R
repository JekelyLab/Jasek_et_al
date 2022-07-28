#code to retrieve the desmosomal connectome from CATMAID and convert it into igraph and tbl_graph format

source("code/Packages_and_Connection.R")

# get info about all desmo connectors in the Naomi project
all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

# extract ids of desmo connectors
all_desmo_ids <- sapply(all_desmo_connectors$connectors, "[[", 1)

# find skids connected to the desmosomes
partners1 <- as.data.frame(sapply(all_desmo_connectors$partners, "[[", 1))
partners1 <- as.vector(sapply(partners1, "[[", 3))

partners2 <- as.data.frame(sapply(all_desmo_connectors$partners, "[[", 2))
partners2 <- as.vector(sapply(partners2, "[[", 3))

#list of all unique skids
all_skids <- unique(c(partners1, partners2))

#generate a graph with all_skids as nodes
g <- make_empty_graph() %>%
  add_vertices(length(all_skids), attr = list(name = all_skids))

#convert to tidy graph
g <- g %>%
  as_tbl_graph()

#add edges based on partner1 and partner2 list with bind_edges, use name as node_key
for (i in 1:length(partners1)){
     g <- g  %>% bind_edges(data.frame(from = as.character(partners1[i]), 
                        to = as.character(partners2[i])), node_key = 'name') }
#list nodes
V(g)$name
#list edges
E(g)

# check connected components  ---------------------------
cl <- components(g)
#the largest subnetwork is membership 1
length(which(cl$membership == 1))

#create a new graph only from the largest connected component
g.largest = induced_subgraph(g, which(cl$membership == 1))

#assign weight of 1 to each edge
E(g.largest)$weight = 1

#sum edges between the same nodes with simplify
desmo_conn_graph <- igraph::simplify(g.largest, remove.loops = FALSE,
            edge.attr.comb = list(weight = "sum", function(x) length(x)) )

E(desmo_conn_graph)

#clustering with Leiden algorithm
partition <- leiden(desmo_conn_graph, weights = E(desmo_conn_graph)$weight, 
                    partition_type = "ModularityVertexPartition",
                    resolution_parameter = 1,
                    n_iterations = -1, seed = 42)
max(partition)

#define colors
blues <- brewer.pal(9, 'Blues')
col=c(brewer.pal(12, 'Paired'), sample(blues, (max(partition)-12), replace=TRUE))

#assign partition value and color to nodes
# If vertices are in the same cluster/partition, give them the same colour
#convert to tidy graph
desmo_conn_graph.tb <- desmo_conn_graph %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(partition = partition) %>%
  mutate(color = col[partition])


# export graph for gephi to do force field layout -------------------------

# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(desmo_conn_graph.tb)
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "data/full_desmo_connectome_graph.gexf")

# read gexf file ----------------------------------------------------------

#read the gephi connectome file with nodes positioned by force field clustering 
#at gexf export, positions were normalised (0,1)

#import gephi file (will be an undirected graph)
conn_gexf <- rgexf::read.gexf("data/full_desmo_connectome_force_layout.gexf")

#get coordinates from imported gephi file
coords <- as_tibble(x = conn_gexf$nodesVizAtt$position$x) %>%
  mutate(x = value) %>%
  mutate(y = conn_gexf$nodesVizAtt$position$y) %>%
  mutate(skid = conn_gexf$nodes$id)
#sort by skid
x_coord <-  as.list(arrange(coords, desc(skid)) %>%
                      select(x) )
y_coord <-  as.list(arrange(coords, desc(skid)) %>%
                      select(y) )

#assign coordinates from imported gephi graph to original CATMAID graph
desmo_conn_graph.tb <- activate(desmo_conn_graph.tb, nodes) %>%
  arrange(desc(name)) %>%  #sort by skid to match coordinate list
  mutate(x = unlist(x_coord),  #add coordinates
         y = unlist(y_coord))


# search for names and annotations --------------------------------------------------

#nodes are named with the skid, use these to get neuron names from CATMAID
#get names with the tidygraph pull function
skids <- pull(desmo_conn_graph.tb, name)
cell_names<-list()
for(i in 1:length(skids)){
  name<-catmaid_get_neuronnames(skids[i], pid = 11)
  cell_names[i] <- name
}

#check duplicated names (there should be none)
cell_names[duplicated(cell_names)]
length(cell_names)

#iterate through skid list, get annotations for the first skid in every cell type
#check for the presence of the annotations and add the annotation to the type_of_cell list
type_of_cell <- list()
annot_to_search <- c("muscle", "epithelia_cell", "basal lamina", 
                     "ciliated cell", "glia cell")
for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      type_of_cell[i] <- annot_to_search[j]
      break()
    } else {type_of_cell[i] <- "other" }
  } 
  print (i)
}

#annotations to search for
annot_to_search <- c("episphere", "segment_0", "segment_1", 
                     "segment_2", "segment_3", "pygidium")
segment_of_cell <- list()
for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      segment_of_cell[i] <- annot_to_search[j]
      break()
    } else {segment_of_cell[i] <- "fragment" }
  } 
  print (i)
}

#annotations to search for
annot_to_search <- c("left_side", "right_side", "middle")
side_of_cell <- list()
for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      side_of_cell[i] <- annot_to_search[j]
      break()
    } else {side_of_cell[i] <- "other" }
  } 
  print (i)
}



#add type of cell to node$group (can be used for visualisation)
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  mutate(group = unlist(type_of_cell)) %>%
  mutate(class = unlist(type_of_cell)) %>%
  mutate(skid = name) %>%
  mutate(name = unlist(cell_names)) %>%
  mutate(segment = unlist(segment_of_cell)) %>%
  mutate(side = unlist(side_of_cell))


# VisNetwork conversion ---------------------------------------------------

{
  #convert to visNetwork graph
  conn_graph.visn <- toVisNetworkData(as.igraph(desmo_conn_graph.tb))  
  
  ## copy column "weight" to new column "value" in list "edges"
  conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
  
  #color will be assigned by group (type of cell)
  #visGroups - $nodes$color takes precedence
  
  #save vis graph with annotations as R data file and txt file printed with dput() to get an exact copy
  saveRDS(conn_graph.visn, "supplements/desmo_connectome_graph.rds")
  writeLines(capture.output(dput(conn_graph.visn)), "supplements/desmo_connectome_graph.txt")
}

# plot graph with coordinates from gephi ----------------------------------

{
  coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol=2)
  
  visNet2 <- visNetwork(conn_graph.visn$nodes,conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=1, max=25),
             color = list(inherit=TRUE, opacity=0.7),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 0.5, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(border='black'),
             opacity = 1, 
             font = list(size = 20)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',labelOnly=FALSE), 
               width = 500, height = 500, autoResize = FALSE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    addFontAwesome() %>%
    visEvents(selectNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'label_long'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].label_long, label_long : label_info[0].label});
            }") %>% 
    visEvents(blurNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'label_long'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].label_long, label_long : label_info[0].label});
  }")
  
  #save as html
  saveNetwork(visNet2, "pictures/Full_desmo_connectome_modules.html", selfcontained = TRUE)
  webshot2::webshot(url="pictures/Full_desmo_connectome_modules.html",
                    file="pictures/Full_desmo_connectome_modules_webshot.png",
                    vwidth = 550, vheight = 550, #define the size of the browser window
                    cliprect = c(50, 60, 480, 480), zoom=8, delay = 5)
}

