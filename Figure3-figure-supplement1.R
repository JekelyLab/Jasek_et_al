#quantification of the modularity of the desmosomal connectome and random scale-free and Erdos graphs of similar statistics
#code for Figure3-figure-supplement1 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely 2021 March

#set working directory

setwd('/workdir/')

#read csv exported from gephi as table
desmo_conn_table <- read.csv('desmosomal-connectome.csv', sep = ";", header = T) 

dim(desmo_conn_table[2:2525])
desmo_conn_matrix <- as.matrix(desmo_conn_table[2:2525],nrow=nrow(desmo_conn_table),ncol=ncol(desmo_conn_table)-1)
desmo_conn_graph <- graph_from_adjacency_matrix(desmo_conn_matrix, mode = "undirected", weighted = T,
                            diag = TRUE, add.colnames = NULL, add.rownames = NA)

###############################
#graph analysis
##############################

# Load packages
library(igraph)
#https://rdrr.io/cran/igraph/man/

#check connected componenets
cl <- components(desmo_conn_graph)
components(desmo_conn_graph)
#create a new graph only from the largest connected component
desmo_conn_graph_largest = induced_subgraph(desmo_conn_graph, which(cl$membership == 1))

#a few iGraph methods
#clustering
modularity(cluster_louvain(desmo_conn_graph))
edge_density(desmo_conn_graph)
#number of edges
length(E(desmo_conn_graph))
#number of nodes
length(V(desmo_conn_graph))
#number of edges
gsize(desmo_conn_graph)
is_weighted(desmo_conn_graph)
sum(E(desmo_conn_graph_largest)$weight)

#convert into binary matrix
desmo_conn_matrix_bi<-as.matrix((desmo_conn_matrix>0)+0)
desmo_conn_graph_bi <- graph_from_adjacency_matrix(desmo_conn_matrix_bi, mode = "undirected", weighted = T,
                                                diag = TRUE, add.colnames = NULL, add.rownames = NA)

#this is how we define an Erdos graph with the same number of nodes and edges as the desmosomal connectome graph
erdos_graph <- erdos.renyi.game(length(V(desmo_conn_graph)), 
          length(E(desmo_conn_graph)),type = "gnm",directed = FALSE,loops = FALSE)

#sample 1000 Erdos graphs with same number of nodes and edges as the desmosomal graph and return their modularity score in a list
modularity_erdos <- lapply(1:1000, function(x) 
  x=modularity(cluster_louvain(erdos.renyi.game(length(V(desmo_conn_graph)), 
  length(E(desmo_conn_graph)),
  type = "gnm",directed = FALSE,loops = FALSE))))

#do the same but assing weights from the desmo graph
modularity_erdos_weighted <- lapply(1:1000, function(x) {
  erdos_graph <- erdos.renyi.game(length(V(desmo_conn_graph)), 
                                  length(E(desmo_conn_graph)),
                                  type = "gnm",directed = FALSE,loops = FALSE)
  E(erdos_graph)$weight <- E(desmo_conn_graph)$weight
  x=modularity(cluster_louvain(erdos_graph))}
  )

#sample subgraphs from the weighted desmosomal graph
modularity_desmo_subgraphs <- lapply(1:1000, function(x)
{vids <- sample(V(desmo_conn_graph),length(V(desmo_conn_graph))-100)
x <- modularity(cluster_louvain(induced_subgraph(desmo_conn_graph, vids, impl = "auto")))}
)

#sample subgraphs from the binarised desmosomal graph
modularity_desmo_bi_subgraphs <- lapply(1:1000, function(x)
{vids <- sample(V(desmo_conn_graph_bi),length(V(desmo_conn_graph))-100)
x <- modularity(cluster_louvain(induced_subgraph(desmo_conn_graph_bi, vids, impl = "auto")))}
)


#generate a scale-free graph

scale_free_graph <- sample_fitness_pl(
  length(V(desmo_conn_graph)),
  length(E(desmo_conn_graph)),
  4,
  exponent.in = -1,
  loops = FALSE,
  multiple = FALSE,
  finite.size.correction = TRUE
)

#sample 1000 scale-free graphs with same number of edges and nodes as the desmosomal graph and return their modularity score in a list
modularity_sf <- lapply(1:1000, function(x) 
  x=modularity(cluster_louvain(scale_free_graph <- sample_fitness_pl(
    length(V(desmo_conn_graph)),
    length(E(desmo_conn_graph)),
    4,
    exponent.in = -1,
    loops = FALSE,
    multiple = FALSE,
    finite.size.correction = TRUE
  ))))

#do the same but assing weights from the desmo graph
modularity_sf_weight <- lapply(1:1000, function(x) {
  sf_graph <- scale_free_graph <- sample_fitness_pl(
    length(V(desmo_conn_graph)),
    length(E(desmo_conn_graph)),
    4,
    exponent.in = -1,
    loops = FALSE,
    multiple = FALSE,
    finite.size.correction = TRUE
  )
  E(sf_graph)$weight <- E(desmo_conn_graph)$weight
  x=modularity(cluster_louvain(sf_graph))
}
)



#plot histograms
pdf(file='Desmosomal_and_Erdos_graph_modularity.pdf', width=10, height = 8)
{
hist_desmo_subgraphs <- hist(as.numeric(modularity_desmo_subgraphs),
                  xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
hist_desmo_bi_subgraphs <- hist(as.numeric(modularity_desmo_bi_subgraphs),
                  xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
hist_erdos <- hist(as.numeric(modularity_erdos),
                  xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
hist_erdos_weighted <- hist(as.numeric(modularity_erdos_weighted),
                  xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
hist_sf <- hist(as.numeric(modularity_sf),
                            xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
hist_sf_weighted <- hist(as.numeric(modularity_sf_weight),
                            xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
plot(hist_erdos,add=F,xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
plot(hist_desmo_subgraphs,add=T)
plot(hist_desmo_bi_subgraphs,add=T)
plot(hist_erdos_weighted,add=T)
plot(hist_sf,add=T)
plot(hist_sf_weighted,add=T)
}
dev.off()


#calculate median modularity values
median(as.numeric(modularity_erdos))
median(as.numeric(modularity_erdos_weighted))
median(as.numeric(modularity_desmo_bi_subgraphs))
median(as.numeric(modularity_desmo_subgraphs))
median(as.numeric(modularity_sf))
median(as.numeric(modularity_sf_weight))
modularity(cluster_louvain(desmo_conn_graph))



#calculate node weights
weights_desmo <- strength(desmo_conn_graph, vids = V(desmo_conn_graph),
                 loops = TRUE)
weights_erdos <- strength(erdos_graph, vids = V(erdos_graph), mode = "all",
                       loops = TRUE)
weights_scale_free=strength(scale_free_graph, vids = V(scale_free_graph), mode = "all",
                       loops = TRUE)

pdf(file='Desmosomal_graph_weight_degree.pdf', width=10, height = 8)
{
hist_desm_degree <- hist(degree(desmo_conn_graph),freq=F,breaks=20,
                         xlab = 'degree',ylab = 'frequency', add=F,main=NA,xlim=c(0,80),ylim=c(0,0.25))
hist_desm_weight <- hist(weights_desmo,freq=F,breaks=100,
                         xlab = 'weights',ylab = 'frequency', add=F,main=NA,xlim=c(0,80),ylim=c(0,0.25))
}
dev.off()


