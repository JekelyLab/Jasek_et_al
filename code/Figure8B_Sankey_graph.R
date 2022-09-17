#R code to generate the Sankey graph in Figure 8 panel B of the 2021 Jasek et al Desmosomal connectome paper

source("code/Packages_and_Connection.R")

MN_names <- c(
  "MNacicX", "MNchae", "vMN1-2",
  "MNbow", "MNwave", "MNbiramous", 
  "MNhose", "MNcrab", "MNladder", 
  "MNspider type1", "MNspider type2", "MNring"
)

annotations_MN <- c(
  "celltype69", "celltype165", "vMN1-2",
  "celltype67", "celltype68", "celltype63", 
  "celltype66", "celltype65", "celltype151",
  "celltype61", "celltype62", "celltype64"
)

annotations_MN_MUSlist = list()
celltypelist = list()

#read cell types
for (i in c(1:12)){
  annotation = paste("annotation:", annotations_MN[i], "$", sep="")
  #read presyn neuron group by annotation
  celltypelist[[i]] <- read.neurons.catmaid(
    paste(
      "annotation:", annotations_MN[i], "$", sep=""
      ), 
    pid=11, 
    fetch.annotations = F
    ) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotations_MN_MUSlist[[i]] <- catmaid_get_annotations_for_skeletons(
    annotation, 
    pid = 11
    )
}

counter <- length(annotations_MN_MUSlist)

#read all non-neuronal celltypes from 37-89 (the muscles) and all annotations
for (i in c(37:89)){
  counter <- counter+1
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation
  celltypelist[[counter]] <- read.neurons.catmaid(
    annotation, 
    pid=11, 
    fetch.annotations = F
    ) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotations_MN_MUSlist[[counter]] <- catmaid_get_annotations_for_skeletons(
    annotation, 
    pid = 11
    )
}

#define empty synapse list with the right dimensions
synapse_list <- vector(
  "list", 
  length(annotations_MN_MUSlist)*length(annotations_MN_MUSlist)
)

#retrieve unique skids for a celltype from the annotation_celltype list
unique(annotations_MN_MUSlist[[1]]$skid)

list_position <- 0
cycle <- 0
for (celltype_skids_pre in annotations_MN_MUSlist){
  presyn_skids <- unique(celltype_skids_pre$skid)
  for (celltype_skids_post in annotations_MN_MUSlist){
    postsyn_skids <- unique(celltype_skids_post$skid)
    assign("celltype_conn", NULL, envir = .GlobalEnv)  #in every iteration we empty the connectivity list  
    # get connectors betwwen neurons of interest
    celltype_conn = catmaid_get_connectors_between(
      pre=presyn_skids, 
      post=postsyn_skids, 
      pid=11
      )
    N_synapses=nrow(celltype_conn)
    if(length(celltype_conn) == 0) {N_synapses <- 0}
    list_position <- list_position + 1
    synapse_list[[list_position]] <- N_synapses
  }
  cycle <- cycle + 1;print (cycle)
}

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(
  unlist(synapse_list), 
  byrow=TRUE, 
  nrow=length(annotations_MN_MUSlist) 
)

#we make a celltype name list using the first cell in every set (not perfect)
celltype_names <- list()
for (df in celltypelist){
  #to retrieves all the neuron names in one element of the celltype list
  neuro_names <- as.character(attr(df,"df")$name)
  celltype_names[[length(celltype_names) + 1]] <- paste(
    neuro_names[1], sep = "_"
    )
}


#truncate the celltype names at the first _
celltype_names <- sub(x=celltype_names, "_(.)*", replacement="")

synapse_matrix <- as.data.frame(synapse_matrix)

#assign column names to matrix
synapse_matrix <- setNames(synapse_matrix, as.character(celltype_names))

#assign row names to matrix
rownames(synapse_matrix) <- as.character(celltype_names)
synapse_matrix <- as.matrix(synapse_matrix)
synapse_matrix

# graph visualization -----------------------------------------------------

celltype_conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"), 
  weighted = TRUE,  
  diag = TRUE, 
  add.colnames = NULL, 
  add.rownames = NA
)

celltype_conn_graph <- celltype_conn_graph %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  filter(weight > 9)

#extract strongly connected componenets
celltype_conn_graph_conn <- components(
  as.igraph(celltype_conn_graph), 
  mode = c("weak")
)

celltype_conn_graph_conn <- induced_subgraph(
  celltype_conn_graph, 
  celltype_conn_graph_conn$members == 1, 
  impl="copy_and_delete"
)

#cluster the graph
members <- leiden(celltype_conn_graph_conn)

# Convert to object suitable for networkD3
celltype_conn_graph_d3 <- igraph_to_networkD3(
  celltype_conn_graph_conn, 
  group = members
)

#The NodeGroup vector in the Nodes data frame needs to be non-numeric so we convert it to character
celltype_conn_graph_d3$nodes$group <- as.character(
  celltype_conn_graph_d3$nodes$group
)


# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain([]) 
.range(["#E69F00", "#56B4E9", "#009E73", 
"#0072B2", "#D55E00", 
"#CC79A7", "#000000"])'

# Plot
Sankey <- sankeyNetwork(Links = celltype_conn_graph_d3$links, 
              Nodes = celltype_conn_graph_d3$nodes, 
              Source = "source",
              Target = "target", 
              NodeID = "name",
              Value = 'value', 
              units = "", 
              NodeGroup="group",
              LinkGroup = NULL,
              colourScale = my_color, 
              fontSize = 15,
              fontFamily = "Arial", 
              nodeWidth = 22, 
              nodePadding = 9, 
              margin = NULL,
              height = 400, 
              width = 700, 
              iterations = 1000, 
              sinksRight = TRUE
)

Sankey

saveNetwork(Sankey, "pictures/Sankey_MN_MUS_circuit.html")
webshot2::webshot(url="pictures/Sankey_MN_MUS_circuit.html",
                  file="pictures/Sankey_MN_MUS_circuit.png",
                  vwidth = 1000, vheight = 600, #define the size of the browser window
                  cliprect = c(0,0,700, 410), zoom = 10, delay = 2)
