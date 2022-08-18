#code to retrieve the desmosomal connectome from CATMAID 
#and convert it into igraph, tbl_graph and visNetwork format
#the first exported gexf graph has to be imported to Gephi for force-field clustering 
#the code also generates a grouped desmosomal graph and saves it in different formats
#final files saved to source data
#these network files are sourced by the coded to generate figures 1, 3 etc.

source("code/Packages_and_Connection.R")

# get info about all desmo connectors in the Naomi project
all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

# extract ids of desmo connectors
all_desmo_ids <- sapply(all_desmo_connectors$connectors, "[[", 1)

#define empty vectors
desmo_id <- c(); x <- c(); y <- c(); z <- c(); partner1 <- c(); partner2 <- c()
#retrieve all desmosomes by id individually and parse their xyz coordinates and partners
for (i in 1:length(all_desmo_ids)) {
  id <- all_desmo_ids[i]
  desmo <- catmaid_fetch(path = paste("11/connectors/", id, sep = ''), 
                         body = list(relation_type="desmosome_with", 
                                     with_partners="true"))
  desmo_id[i] <- desmo$connector_id
  x[i] <- desmo$x
  y[i] <- desmo$y
  z[i] <- desmo$z
  partner1[i] <- desmo$partners[[1]]$skeleton_id
  partner2[i] <- desmo$partners[[2]]$skeleton_id
}

#assemble all in a tibble
desmo_with_partners <- tibble('desmo_id' = desmo_id,
                              'x' = x,
                              'y' = y,
                              'z' = z,
                              'partner1' = partner1,
                              'partner2' = partner2
)


desmo_with_partners

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
annot_to_search <- c("notopodium", "neuropodium")
parapod_region <- list()
for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      parapod_region[i] <- annot_to_search[j]
      break()
    } else {parapod_region[i] <- "non_parapodial" }
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
  mutate(side = unlist(side_of_cell)) %>%
  mutate(parapod_region = unlist(parapod_region))


#search celltype_non_neuronal annotations
celltype_id <- list()
celltype_id_num <- list()
annot_to_search <- c("fragmentum", "basal lamina", 
                     paste("celltype_non_neuronal", c(1:21), sep = ""),
                     "chaeta_neuro", "chaeta_noto",
                     "acicula_neuro", "acicula_noto",
                     "circumacicular_neuro", "circumacicular_noto",
                     "hemichaetal_neuro", "hemichaetal_noto",
                     "ER-circumchaetal_neuro", "ER-circumchaetal_noto",
                     "noER-circumchaetal_neuro", "noER-circumchaetal_noto",
                     "circumchaetal_neuro", "circumchaetal_noto",
                     paste("celltype_non_neuronal", c(28:92), sep = ""),
                     "with_soma")

#detailed annotations to search for, separate for notopodium and neuropodium for the acicular/chaetal complex
#celltype_non_neuronal22 chaeta chaeta_neuro chaeta_noto
#celltype_non_neuronal23 acicula acicula_neuro acicula_noto
#celltype_non_neuronal24 circumacicular circumacicular_neuro circumacicular_noto
#celltype_non_neuronal25 hemichaetal hemichaetal_neuro hemichaetal_noto
#celltype_non_neuronal26 ER-circumchaetal ER-circumchaetal_neuro ER-circumchaetal_noto
#celltype_non_neuronal27 noER-circumchaetal noER-circumchaetal_neuro noER-circumchaetal_noto
#cellgroup_non-neuronal1_circumchaetal circumchaetal_neuro circumchaetal_noto

#search annotations and assign to celltype_id list (annotations will be
#progressively overwritten if duplicate)

for(i in seq_along(skids)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      celltype_id[i] <- annot_to_search[j]
      celltype_id_num[i] <- j
      break()
    } else {celltype_id[i] <- "other"
            celltype_id_num[i] <- length(annot_to_search)+1 }
  } 
  print (i)
}

#add type of cell to node$group (can be used for visualisation)
desmo_conn_graph.tb <- desmo_conn_graph.tb %>%
  mutate(celltype = unlist(celltype_id)) %>%
  mutate(celltype_num = unlist(celltype_id_num))

celltype_id_list[celltype_id_list == "other"]

#add node weighted degree
desmo_conn_graph <- as.igraph(desmo_conn_graph.tb)
V(desmo_conn_graph)$weighted_degree <- strength(desmo_conn_graph, vids = V(desmo_conn_graph),
                                       mode = 'all', loop = TRUE)
desmo_conn_graph.tb <- as_tbl_graph(desmo_conn_graph)

#save in igraph format
saveRDS(as.igraph(desmo_conn_graph.tb), "source_data/Figure3_source_data1.rds")
#read the saved igraph format graph file from supplements/
desmo_conn_graph <- readRDS("source_data/Figure3_source_data1.rds")


# contract vertices by cell type to make grouped graph ---------------------------------
mapping_df <- data.frame(as_tibble(desmo_conn_graph.tb %>%
            select(celltype_num) ))[,1]

desmo_grouped_graph <- contract.vertices(desmo_conn_graph, 
                                         mapping = mapping_df,
                  vertex.attr.comb = list(name = "first", 
                                          "first"))
#sum edges between the same nodes with simplify
desmo_grouped_graph <- igraph::simplify(desmo_grouped_graph, 
                          remove.loops = FALSE,
                          edge.attr.comb = list(degree = "sum", 
                                  function(x) length(x)) )

desmo_grouped_graph.tb <- desmo_grouped_graph %>%
  as_tbl_graph() %>%
  select(celltype, CATMAID_name, class, weighted_degree) %>%
  filter(class != "basal lamina") %>%
  filter(class != "other") %>%
#  filter(class != "epithelia_cell") %>%
  filter(CATMAID_name != "character(0)") %>%
  filter(CATMAID_name != "with_soma")

#get the largest connected component
desmo_grouped_graph <- as.igraph(desmo_grouped_graph.tb)
V(desmo_grouped_graph)$class
cl <- components(desmo_grouped_graph)
#the largest subnetwork is membership 1
length(which(cl$membership == 1))

#create a new graph only from the largest connected component
desmo_grouped_graph = induced_subgraph(desmo_grouped_graph, which(cl$membership == 1))
desmo_grouped_graph <- desmo_grouped_graph %>% 
  as_tbl_graph() %>%
  mutate(celltype = unlist(celltype)) %>%
  mutate(CATMAID_name = unlist(CATMAID_name)) %>%
  mutate(class = unlist(class)) %>%
  mutate(CATMAID_name = sub("_.*$", "", CATMAID_name)) %>% #truncate names at _ and \s
  mutate(CATMAID_name = sub("\\s.*$", "", CATMAID_name))

E(desmo_grouped_graph)$weight
desmo_grouped_graph %>%
  as_tbl_graph()

#add noto/neuro to some node names
desmo_grouped_graph <- desmo_grouped_graph %>%
  as_tbl_graph() %>%
  mutate(CATMAID_name = ifelse(celltype == "chaeta_neuro", 
                               "chaeta_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "chaeta_noto", 
                               "chaeta_noto",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "acicula_neuro", 
                               "acicula_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "acicula_noto", 
                               "acicula_noto",
                               CATMAID_name))  %>%
  mutate(CATMAID_name = ifelse(celltype == "circumacicular_neuro", 
                               "circumacicular_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "circumacicular_noto", 
                               "circumacicular_noto",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "circumchaetal_neuro", 
                               "circumchaetal_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "circumchaetal_noto", 
                               "circumchaetal_noto",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "ER-circumchaetal_neuro", 
                               "ER-circumchaetal_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "ER-circumchaetal_noto", 
                               "ER-circumchaetal_noto",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "noER-circumchaetal_neuro", 
                               "noER-circumchaetal_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "noER-circumchaetal_noto", 
                               "noER-circumchaetal_noto",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "hemichaetal_neuro", 
                               "hemichaetal_neuro",
                               CATMAID_name)) %>%
  mutate(CATMAID_name = ifelse(celltype == "hemichaetal_noto", 
                               "hemichaetal_noto",
                               CATMAID_name))

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

desmo_grouped_graph
saveRDS(desmo_grouped_graph, "source_data/Figure3_figure_supplement1_source_data3.rds")
writeLines(capture.output(dput(desmo_grouped_graph)), "source_data/Figure3_figure_supplement1_source_data3.txt")

}

