#R code to generate the Sankey graph in Figure 8 panel B of the 2021 Jasek et al Desmosomal connectome paper

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#load RColorBrewer for color palette 
library(RColorBrewer)

#for regex
library(dplyr)
library(readr)

#set working directory
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Muscles/Figures/Figure8-desmo-MN/')
#setwd("/working_directory/")

# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info


MN_names <- c("MNacicX","MNchae","vMN1-2","MNbow", "MNwave","MNbiramous", "MNhose", "MNcrab", "MNladder", "MNspider type1", "MNspider type2","MNring")
annotations_MN <- c("celltype69","celltype165","vMN1-2","celltype67", "celltype68", "celltype63", "celltype66", "celltype65", "celltype151","celltype61","celltype62","celltype64")
length(MN_names)

annotations_MN_MUSlist = list()
celltypelist = list()

#first we read all celltypes from 1-182 and all annotations
for (i in c(1:12)){
  annotation = paste("annotation:", annotations_MN[i], "$", sep="")
  #read presyn neuron group by annotation
  celltypelist[[i]] <- read.neurons.catmaid(paste("annotation:", annotations_MN[i], "$", sep=""), pid=11, fetch.annotations = F) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotations_MN_MUSlist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

counter=length(annotations_MN_MUSlist)
#read all non-neuronal celltypes from 37-89 (the muscles) and all annotations
for (i in c(37:89)){
  counter <- counter+1
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation
  celltypelist[[counter]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = F) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotations_MN_MUSlist[[counter]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define empty synapse list with the right dimensions
synapse_list <- vector("list", length(annotations_MN_MUSlist)*length(annotations_MN_MUSlist))
#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids

#retrieve unique skids for a celltype from the annotation_celltype list
unique(annotations_MN_MUSlist[[1]]$skid)

list_position=0;cycle=0
for (celltype_skids_pre in annotations_MN_MUSlist){
  presyn_skids <- unique(celltype_skids_pre$skid)
  for (celltype_skids_post in annotations_MN_MUSlist){
    postsyn_skids <- unique(celltype_skids_post$skid)
    assign("celltype_conn", NULL, envir = .GlobalEnv)  #in every iteration we empty the connectivity list  
    # get connectors betwwen neurons of interest
    celltype_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
    N_synapses=nrow(celltype_conn)
    if(length(celltype_conn) == 0) {N_synapses <- 0}
    list_position <- list_position + 1
    synapse_list[[list_position]] <- N_synapses
  }
  cycle <- cycle + 1;print (cycle)
}

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=length(annotations_MN_MUSlist) )

#we make a celltype name list using the first cell in every set (not perfect)
celltype_names=list()
for (df in celltypelist){
  #to retrieves all the neuron names in one element of the celltype list
  neuro_names <- as.character(attr(df,"df")$name)
  celltype_names[[length(celltype_names) + 1]] <- paste(neuro_names[1], sep = "_")
}


#truncate the celltype names at the first _
celltype_names <- sub(x=celltype_names, "_(.)*", replacement="")

synapse_matrix = as.data.frame(synapse_matrix)

#assign column names to matrix
synapse_matrix=setNames(synapse_matrix, as.character(celltype_names))

#assign row names to matrix
rownames(synapse_matrix) <- as.character(celltype_names)
synapse_matrix = as.matrix(synapse_matrix)



###############################
#graph visualisation
##############################

# Load igraph
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(networkD3)
library(htmlwidgets)

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
celltype_conn_graph <- graph_from_adjacency_matrix(synapse_matrix,
                                                   mode = c("directed"),
                                                   weighted = NULL,  diag = TRUE, add.colnames = NULL, add.rownames = NA)

#extract strongly connected componenets
celltype_conn_graph_conn <- components(celltype_conn_graph, mode = c("strong"))
celltype_conn_graph_conn <- induced_subgraph(celltype_conn_graph, celltype_conn_graph_conn$members, impl="copy_and_delete")

#cluster the graph
wc <- cluster_walktrap(celltype_conn_graph_conn)
members <- membership(wc)

# Convert to object suitable for networkD3
celltype_conn_graph_d3 <- igraph_to_networkD3(celltype_conn_graph_conn, group = members)

#The NodeGroup vector in the Nodes data frame needs to be non-numeric so we convert it to character
celltype_conn_graph_d3$nodes$group <- as.character(celltype_conn_graph_d3$nodes$group)

# Create force directed network plot
forceNetwork(Links = celltype_conn_graph_d3$links, Nodes = celltype_conn_graph_d3$nodes, 
             Source = 'source', Target = 'target', 
             NodeID = 'name', Group = 'group')

#we need to define the value for proper plotting
celltype_conn_graph_d3$links$value=1

#we can reassign the group if we want to color by MUS vs MN instead of network group (optional)
celltype_conn_graph_d3$nodes$group=c(c(rep('1',times=12)),c(rep('2',times=42)))

# Plot
sankeyNetwork(Links = celltype_conn_graph_d3$links, Nodes = celltype_conn_graph_d3$nodes, Source = "source",
              Target = "target", NodeID = "name",Value = 'value', 
              units = "", NodeGroup="group",
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"), fontSize = 10,
              fontFamily = NULL, nodeWidth = 40, nodePadding = 4, margin = NULL,
              height = 500, width = 600, iterations = 32, sinksRight = TRUE)


#save it by first opening it 'as a webpage a'in a new window' and then save from browser as pdf = pdf() saving from R did not work
