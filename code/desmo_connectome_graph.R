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

E(desmo_conn_graph)$weight

#clustering with Leiden algorithm
partition <- leiden(desmo_conn_graph, weights = E(desmo_conn_graph)$weight, 
                    partition_type = "RBConfigurationVertexPartition",
                    resolution_parameter = 0.27,
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

#Force field clustering was carried out with the Force Atlas tool in Gephi (0.9.2)
#The inertia was set to 0.1, repulsion strength was 35, 
#attraction strength was 10, maximum displacement was 5, 
#gravity was 50, speed was 5 and the 'attraction distribution' option was selected. 
#The 'auto stabilise function' was off. Towards the end of the clustering 
#the 'adjust by sizes' option was also selected. 
#To prevent node overlap, we then run the 'Noverlap' function. 

# read gexf file ----------------------------------------------------------

#read the gephi connectome file with nodes positioned by force field clustering 
#at gexf export, positions were normalised (0,1)

#import gephi file (will be an undirected graph)
conn_gexf <- rgexf::read.gexf("data/full_desmo_connectome_force_layout.gexf")

#get coordinates from imported gephi file
coords <- as_tibble(x = conn_gexf$nodesVizAtt$position$x) %>%
  rename('x' = value) %>%
  mutate(y = conn_gexf$nodesVizAtt$position$y) %>%
  mutate(skid = as.numeric(conn_gexf$nodes$id))
coords

#sort by skid
x_coord <-  as.list(arrange(coords, desc(skid)) %>%
                      select(x) )
y_coord <-  as.list(arrange(coords, desc(skid)) %>%
                      select(y) )

#assign coordinates from imported gephi graph to original CATMAID graph
desmo_conn_graph.tb <- activate(desmo_conn_graph.tb, nodes) %>%
  arrange(desc(as.numeric(name))) %>%  #sort by skid to match coordinate list
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
                     "acicula", "chaeta", "circumacicular",
                     "circumchaetal", "ciliated cell",
                     "hemichaetal")
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
    } else {segment_of_cell[i] <- "fragmentum" }
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
    } else {side_of_cell[i] <- "fragmentum" }
  } 
  print (i)
}

#add type of cell to node$group (can be used for visualisation)
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  mutate(group = unlist(type_of_cell)) %>%
  mutate(class = unlist(type_of_cell)) %>%
  mutate(CATMAID_name = unlist(cell_names)) %>%
  mutate(segment = unlist(segment_of_cell)) %>%
  mutate(side = unlist(side_of_cell))

#search celltype_non_neuronal annotations
celltype_id <- list()
annot_to_search <- c("fragmentum", "basal lamina", 
                     paste("celltype_non_neuronal", c(1:92), sep = ""),
                     "with_soma")

#search annotations and assign to celltype_id list (annotations will be
#progressively overwritten if duplicate)

for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      celltype_id[i] <- annot_to_search[j]
      break()
    } else {celltype_id[i] <- "other" }
  } 
  print (i)
}

#add type of cell to node$group (can be used for visualisation)
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  mutate(celltype = unlist(celltype_id))
#if there are still nodes with onle 'other' annotation, delete these
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  filter(celltype != 'other')

#check unique celltype ids
N_unique_celltype_ids <- dim(unique(as_tibble(desmo_conn_graph.tb %>%
                                                select(celltype))) 
)[1]

#make a list and number celltype ids (will need for vertex contraction)
celltype_id_list <- as.data.frame(unique(as_tibble(desmo_conn_graph.tb %>%
                                                     select(celltype))) %>%
                                    mutate(grouped_id = c(1:N_unique_celltype_ids)) )
#run again cell type annotation search and add id number as well
#search celltype_non_neuronal annotations
celltype_id <- list()
celltype_id_num <- list()
annot_to_search <- rev(celltype_id_list[,1])
#check if 'other' is in the annotation list (should not be)
celltype_id_list[celltype_id_list == "other"]

for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      celltype_id[i] <- celltype_id_list[j,1]
      celltype_id_num[i] <- celltype_id_list[j,2]
      break()
    } else {celltype_id[i] <- "other"
    celltype_id_num[i] <- length(celltype_id_list[,1])+1 }
  } 
  print (i)
}

#add type of cell to node$group (can be used for visualisation)
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  mutate(celltype = unlist(celltype_id)) %>%
  mutate(celltype_num = unlist(celltype_id_num))

#save in igraph format
saveRDS(as.igraph(desmo_conn_graph.tb), "source_data/Figure3_source_data1.rds")
#read the saved igraph format graph file from supplements/
desmo_conn_graph <- readRDS("source_data/Figure3_source_data1.rds")


# contract vertices by cell type to make grouped graph ---------------------------------

mapping_df <- data.frame(as_tibble(desmo_conn_graph.tb %>%
            select(celltype_num) ))[,1]

V(desmo_conn_graph)$weight <- strength(desmo_conn_graph, 
                                       vids = V(desmo_conn_graph),
                                       mode = 'all',
                                       loops = TRUE)

desmo_grouped_graph <- contract.vertices(desmo_conn_graph, 
                                         mapping = mapping_df,
                  vertex.attr.comb = list(x="max", 
                                          y = "max",
                                          name = "first", 
                                          "first"))

#sum edges between the same nodes with simplify
desmo_grouped_graph <- igraph::simplify(desmo_grouped_graph, 
                          remove.loops = FALSE,
                          edge.attr.comb = list(weight = "sum", 
                                  function(x) length(x)) )
V(desmo_grouped_graph)$class
desmo_grouped_graph.tb <- desmo_grouped_graph %>%
  as_tbl_graph() %>%
  select(celltype, CATMAID_name, class, weight) %>%
  filter(class != "basal lamina") %>%
  filter(class != "other") %>%
#  filter(class != "epithelia_cell") %>%
  filter(CATMAID_name != "character(0)")

#get the largest connected component
desmo_grouped_graph <- as.igraph(desmo_grouped_graph.tb)
V(desmo_grouped_graph)$class
cl <- components(desmo_grouped_graph)
#the largest subnetwork is membership 1
length(which(cl$membership == 1))

#create a new graph only from the largest connected component
desmo_grouped_graph = induced_subgraph(desmo_grouped_graph, which(cl$membership == 1))

# VisNetwork conversion ---------------------------------------------------

{
#convert to visNetwork graph
conn_grouped_graph.visn <- toVisNetworkData(desmo_grouped_graph)  
  
## copy column "weight" to new column "value" in list "edges"
conn_grouped_graph.visn$edges$value <- conn_grouped_graph.visn$edges$weight
conn_graph.visn <- toVisNetworkData(desmo_conn_graph)
#save both ungrouped and grouped vis graph with annotations as R data file and txt file printed with dput() to get an exact copy
saveRDS(conn_graph.visn, "source_data/Figure3_source_data2.rds")
writeLines(capture.output(dput(conn_graph.visn)), "source_data/Figure3_source_data2.txt")

saveRDS(conn_grouped_graph.visn, "source_data/Figure3_figure_supplement1_source_data1.rds")
writeLines(capture.output(dput(conn_grouped_graph.visn)), "source_data/Figure3_figure_supplement1_source_data1.txt")

#read the saved visNetwork file from supplements/
conn_grouped_graph.visn <- readRDS("source_data/Figure3_figure_supplement1_source_data1.rds")

}

# plot grouped graph ----------------------------------

#overwrite group value (partition) with side value or other value (for colouring)
conn_grouped_graph.visn$nodes$group <-  as.character(conn_grouped_graph.visn$nodes$class)

## copy column "weight" to new column "value" in list "edges"
conn_grouped_graph.visn$edges$value <- sqrt(conn_grouped_graph.visn$edges$weight)
## copy column "weight" to new column "value" in list "nodes"
conn_grouped_graph.visn$nodes$value <- conn_grouped_graph.visn$nodes$weight
#for plotting with different colors, remove colour info (which takes precedence over group colour)
conn_grouped_graph.visn$nodes$color <- c()
Reds <- brewer.pal(9, 'Reds')

{
  
#shorten names after _ and \s
conn_grouped_graph.visn$nodes$CATMAID_name  <- sub("_.*$", "", conn_grouped_graph.visn$nodes$CATMAID_name)
conn_grouped_graph.visn$nodes$CATMAID_name  <- sub("\\s.*$", "", conn_grouped_graph.visn$nodes$CATMAID_name)

#assign name to label
conn_grouped_graph.visn$nodes$label <-  conn_grouped_graph.visn$nodes$CATMAID_name
conn_grouped_graph.visn$nodes
visNet <- visNetwork(conn_grouped_graph.visn$nodes, 
                     conn_grouped_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE) %>%
  visPhysics(solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -150000,
                                       centralGravity = 10)) %>%
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
            opacity=0.5) %>%
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
                                     opacity=0.8) 

visNet
  
#save as html
saveNetwork(visNet, "pictures/Grouped_desmo_connectome.html", selfcontained = TRUE)
#create png snapshot
webshot2::webshot(url="pictures/Grouped_desmo_connectome.html",
                    file="pictures/Grouped_desmo_connectome.png",
                    vwidth = 1500, vheight = 1500, #define the size of the browser window
                    cliprect = c(50, 60, 1500, 1500), zoom=5, delay = 2)
}


