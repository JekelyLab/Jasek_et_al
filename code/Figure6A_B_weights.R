#quantification and visualisation of node degrees in the desmosomal connectome 
#code for Figure3-figure-supplement1 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely 2021 March

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# Load packages
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(natverse)
#all methods:
#https://rdrr.io/github/natverse/rcatmaid/man/

# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
source("~/R/conn.R")


#read csv exported from gephi as table
#this is also in the github page as desmosomal-connectome.csv
desmo_conn_table <- read.csv('data/desmosomal-connectome.csv', sep = ";", header = T) 

#read csv exported from gephi containing skids and neuron names
data_lab_table <- read.csv('data/Full_desmosomal_graph_dataLlaboratory.csv', sep = ";", header = T) 

skids <- unlist(data_lab_table[2])
length(skids)

dim(desmo_conn_table[2:2525])
desmo_conn_matrix <- as.matrix(desmo_conn_table[2:2525],nrow=nrow(desmo_conn_table),ncol=ncol(desmo_conn_table)-1)
dim(desmo_conn_matrix)

colnames(desmo_conn_matrix) <- skids
rownames(desmo_conn_matrix) <- skids

desmo_conn_graph <- graph_from_adjacency_matrix(desmo_conn_matrix, mode = "undirected", weighted = T,
                                                diag = TRUE, add.colnames = NULL, add.rownames = NA)


###############################
#graph analysis
##############################

edge_density(desmo_conn_graph)
#number of edges
length(E(desmo_conn_graph))
#number of nodes
length(V(desmo_conn_graph))
#number of edges
gsize(desmo_conn_graph)
#number of edges
is_weighted(desmo_conn_graph)

node_weights <- strength(desmo_conn_graph, weights = NULL)

node_weights_sorted <- sort(node_weights, decreasing = TRUE)

max(node_weights_sorted)

conn <- source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))

high_weight_cells <- nlapply(read.neurons.catmaid(names(node_weights_sorted[1:1000]), pid=11), 
                             function(x) smooth_neuron(x, sigma=6000))


######################
#Load catmaid reference volumes and cells
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                              invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)
acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

##########################################
# 3d plotting 
#########################################

nopen3d() # opens apannable 3d window
mfrow3d(1, 1)  #defines the two scenes
par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
par3d(zoom=0.52)

clear3d()


##############################
#plot cells coloured by weighted degree

for (i in 1:1000){
  alpha_factor=(round(node_weights_sorted[i]/8.3))/10
   print (alpha_factor, " ", i)
  plot3d(high_weight_cells[i], WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=alpha_factor,
         col=hcl.colors(1, palette='Reds 3'))
}

plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 

plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 

plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 

#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_highest_weight_ventral.png")

#lateral view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.55)
rgl.snapshot("pictures/desmosomal_connectome_highest_weight_right.png")

