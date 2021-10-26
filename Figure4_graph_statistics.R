#R code to quantify modularity of the desmosomal connectome, code for Figure 4 of the 2021 Jasek et al desmosomal connectome paper
#quantification of the modularity of the desmosomal connectome and random scale-free and Erdos graphs of similar statistics
#Gaspar Jekely 2021 March


# Load packages
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(leiden)
#https://github.com/TomKellyGenetics/leiden

library(reticulate)

#set working directory

setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Muscles/Figures/Figure4-network-statistics/')


#read csv exported from gephi as table
#this is also in the github page as desmosomal-connectome.csv
desmo_conn_table <- read.csv('desmosomal-connectome.csv', sep = ";", header = T) 

dim(desmo_conn_table[2:2525])
desmo_conn_matrix <- as.matrix(desmo_conn_table[2:2525],nrow=nrow(desmo_conn_table),ncol=ncol(desmo_conn_table)-1)
desmo_conn_graph <- graph_from_adjacency_matrix(desmo_conn_matrix, mode = "undirected", weighted = T,
                            diag = TRUE, add.colnames = NULL, add.rownames = NA)


###############################
#graph analysis
##############################

#check connected components
cl <- components(desmo_conn_graph)
components(desmo_conn_graph)

#clustering with Leiden algorith, run with ModularityVertexPartition and resolution parameter
partition <- leiden(desmo_conn_graph, partition_type = "ModularityVertexPartition", resolution_parameter = 0.95)
#calculate modularity
modularity(desmo_conn_graph,partition)

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
#number of edges
is_weighted(desmo_conn_graph)
sum(E(desmo_conn_graph_largest)$weight)

#convert into binary matrix
desmo_conn_matrix_bi<-as.matrix((desmo_conn_matrix>0)+0)
desmo_conn_graph_bi <- graph_from_adjacency_matrix(desmo_conn_matrix_bi, mode = "undirected", weighted = T,
                                                diag = TRUE, add.colnames = NULL, add.rownames = NA)

#this is how we define an Erdos graph with the same number of nodes and edges as the desmosomal connectome graph
erdos_graph <- erdos.renyi.game(length(V(desmo_conn_graph)), 
          length(E(desmo_conn_graph)),type = "gnm",directed = FALSE,loops = FALSE)

partition <- leiden(desmo_conn_graph_largest, partition_type = "RBConfigurationVertexPartition",
                                            resolution_parameter = 1)
max(partition)
#calculate modularity

system.time(modularity(desmo_conn_graph,leiden(desmo_conn_graph, resolution_parameter = 0.95)))
system.time(modularity(cluster_louvain(desmo_conn_graph)))
#will use louvain clustering to calculate modularity


#define a modularity_leiden function
modularity_leiden <- function(graph) {
  modularity(graph,leiden(graph, resolution_parameter = 0.95))
}


#sample 1000 Erdos graphs with same number of nodes and edges as the desmosomal graph and return their modularity score in a list
modularity_erdos <-lapply(1:10, function(x) 
  modularity_leiden(erdos.renyi.game(length(V(desmo_conn_graph)), length(E(desmo_conn_graph)),type = "gnm",directed = FALSE,loops = FALSE))
)

#do the same but assign weights from the desmo graph
modularity_erdos_weighted <- lapply(1:10, function(x) {
  erdos_graph <- erdos.renyi.game(length(V(desmo_conn_graph)), 
                                  length(E(desmo_conn_graph)),
                                  type = "gnm",directed = FALSE,loops = FALSE)
  E(erdos_graph)$weight <- E(desmo_conn_graph)$weight
  x=modularity_leiden(erdos_graph)}
  )

#calculate the global transitivity (clustering coefficient) of 1000 Erdos graphs
transitivity_erdos <- lapply(1:10, function(x) 
  x=transitivity(erdos.renyi.game(length(V(desmo_conn_graph)), 
  length(E(desmo_conn_graph)),
  type = "gnm",directed = FALSE,loops = FALSE)
))

#sample subgraphs from the weighted desmosomal graph
modularity_desmo_subgraphs <- lapply(1:10, function(x)
{vids <- sample(V(desmo_conn_graph),length(V(desmo_conn_graph))-100)
x <- modularity_leiden(induced_subgraph(desmo_conn_graph, vids, impl = "auto"))}
)

#sample subgraphs from the binarised desmosomal graph
modularity_desmo_bi_subgraphs <- lapply(1:10, function(x)
{vids <- sample(V(desmo_conn_graph_bi),length(V(desmo_conn_graph))-100)
x <- modularity_leiden(induced_subgraph(desmo_conn_graph_bi, vids, impl = "auto"))}
)

#calculate the global transitivity (clustering coefficient) of 1000 subsamples binarised desmosomal graphs 
transitivity_desmo <- lapply(1:10, function(x) 
{vids <- sample(V(desmo_conn_graph_bi),length(V(desmo_conn_graph))-100)
x <- transitivity(induced_subgraph(desmo_conn_graph_bi, vids, impl = "auto"))}
)

#######################################
#generate and quantify scale-free graphs

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
modularity_sf <- lapply(1:10, function(x) 
  x=modularity_leiden(sample_fitness_pl(
    length(V(desmo_conn_graph)),
    length(E(desmo_conn_graph)),
    4,
    exponent.in = -1,
    loops = FALSE,
    multiple = FALSE,
    finite.size.correction = TRUE
)))

#calculate the global transitivity (clustering coefficient) of 1000 scale free graphs
transitivity_sf <- lapply(1:10, function(x) x=transitivity(sample_fitness_pl(
  length(V(desmo_conn_graph)),length(E(desmo_conn_graph)),
  4,exponent.in = -1, loops = FALSE,multiple = FALSE,
  finite.size.correction = TRUE
)))

#do the same but assign weights from the desmo graph
modularity_sf_weight <- lapply(1:10, function(x) {
  sf_graph <- sample_fitness_pl(
    length(V(desmo_conn_graph)),
    length(E(desmo_conn_graph)),
    4,
    exponent.in = -1,
    loops = FALSE,
    multiple = FALSE,
    finite.size.correction = TRUE
  )
  E(sf_graph)$weight <- E(desmo_conn_graph)$weight
  x=modularity_leiden(sf_graph)
}
)

#calculate the length of maximal cliques of 1000 preferential attachment graphs
#run with min=3, was also run with min=4 which returned 0 clique
max_cliques_sf <- lapply(1:10, function(x) x=length(max_cliques(sample_fitness_pl(
  length(V(desmo_conn_graph)),length(E(desmo_conn_graph)),
  4,exponent.in = -1, loops = FALSE,multiple = FALSE,
  finite.size.correction = TRUE
), min=3))
)
sum(as.numeric(max_cliques_sf))
median(as.numeric(max_cliques_sf))
hist(as.numeric(max_cliques_sf))



#####################################
#graphs with preferential attachment

graph_pa_age <- sample_pa_age(
  length(V(desmo_conn_graph)),
  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,
  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,
  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL
)
modularity_leiden(graph_pa_age)
transitivity(graph_pa_age)
max(degree(graph_pa_age))
gsize(graph_pa_age)
gsize(desmo_conn_graph)
length(V(desmo_conn_graph))

#sample 1000 preferential attachment graphs with same number of nodes and very similar number of edges as the desmosomal graph and return their modularity score in a list
modularity_pa <- lapply(1:10, function(x) 
  x=modularity_leiden(graph_pa_age <- sample_pa_age(
    length(V(desmo_conn_graph)),
    pa.exp=1,aging.exp=-2,
    m = 2, aging.bin = 300,
    out.dist = NULL,out.seq = NULL,
    out.pref = FALSE,directed = F,
    zero.deg.appeal = 1,zero.age.appeal = 0,
    deg.coef = 1,age.coef = 1,time.window = NULL
  ))
)

#calculate the global transitivity (clustering coefficient) of 1000 preferential attachment graphs
transitivity_pa <- lapply(1:10, function(x) x=transitivity(graph_pa_age <- sample_pa_age(
  length(V(desmo_conn_graph)),
  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,
  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,
  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL
)))

#calculate the length of maximal cliques of 1000 preferential attachment graphs
#run with min=3, was also run with min=4 which returns no cliques
max_cliques_pa <- lapply(1:10, function(x) x=length(max_cliques(sample_pa_age(
  length(V(desmo_conn_graph)),
  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,
  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,
  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL), min=3))
)
max(as.numeric(max_cliques_pa))
median(as.numeric(max_cliques_pa))
hist(as.numeric(max_cliques_pa))

#plot histograms of modularity scores
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
hist_pa <- hist(as.numeric(modularity_pa),
                         xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
plot(hist_erdos,add=F,xlim=c(0.5,0.9),ylim=c(0,300),ylab = 'count',xlab='modularity', main=NA)
plot(hist_desmo_subgraphs,add=T)
plot(hist_desmo_bi_subgraphs,add=T)
plot(hist_erdos_weighted,add=T)
plot(hist_sf,add=T)
plot(hist_sf_weighted,add=T)
plot(hist_pa,add=T)
}
dev.off()


#calculate median modularity values
median(as.numeric(modularity_erdos))
median(as.numeric(modularity_erdos_weighted))
median(as.numeric(modularity_desmo_bi_subgraphs))
median(as.numeric(modularity_desmo_subgraphs))
median(as.numeric(modularity_sf))
median(as.numeric(modularity_sf_weight))
median(as.numeric(modularity_pa))
modularity_leiden(desmo_conn_graph)

#check package versions
packageVersion("igraph") 
packageVersion("leiden") 


#############################################
#plot histograms of transitivity scores
pdf(file='Desmosomal_Erdos_and_sf_graph_transitivity.pdf', width=10, height = 8)
{
  hist_desmo_transitivity <- hist(as.numeric(transitivity_desmo),
                               xlim=c(0,0.1),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA)
  hist_erdos_transitivity <- hist(as.numeric(transitivity_erdos),
                                  xlim=c(0,0.1),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA,col='red')
  hist_sf_transitivity <- hist(as.numeric(transitivity_sf),
                     xlim=c(0,0.1),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA, col='cyan')
  hist_pa_transitivity <- hist(as.numeric(transitivity_pa),
                               xlim=c(0,0.1),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA, col='cyan')
  
  plot(hist_desmo_transitivity,add=F,xlim=c(0,0.1),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA)
  plot(hist_erdos_transitivity,add=T, col='red')
  plot(hist_sf_transitivity,add=T, col='cyan')
  plot(hist_pa_transitivity,add=T, col='grey50')
 }
dev.off()

median(as.numeric(transitivity_desmo))
median(as.numeric(transitivity_erdos))
median(as.numeric(transitivity_sf))
median(as.numeric(transitivity_pa))

#calculate node weights
weights_desmo <- strength(desmo_conn_graph, vids = V(desmo_conn_graph),
                 loops = TRUE)
weights_erdos <- strength(erdos_graph, vids = V(erdos_graph), mode = "all",
                       loops = TRUE)
weights_scale_free=strength(scale_free_graph, vids = V(scale_free_graph), mode = "all",
                       loops = TRUE)

#other network parameters
mean(weights_desmo)
diameter(desmo_conn_graph_largest)
diameter(erdos_graph)
diameter(scale_free_graph)
diameter(graph_pa_age)

mean_distance(desmo_conn_graph)
mean_distance(erdos_graph)
mean_distance(scale_free_graph)
length(which((coreness(desmo_conn_graph,mode='all'))==5))
max(weights_desmo)
sort(degree(desmo_conn_graph), decreasing=T)
count_motifs(desmo_conn_graph, 4)
count_motifs(erdos_graph, 4)
count_motifs(scale_free_graph, 4)
max(count_triangles(desmo_conn_graph))
(triangles(desmo_conn_graph))

#triangles
count_triangles(desmo_conn_graph)
sort(count_triangles(desmo_conn_graph), decreasing=T)
max(count_triangles(desmo_conn_graph))
max(count_triangles(erdos_graph))
max(count_triangles(scale_free_graph))
max(count_triangles(graph_pa_age))
count_triangles(desmo_conn_graph, vids=V(desmo_conn_graph)[223])
V(desmo_conn_graph)[223]
count_triangles(erdos_graph)

#cliques
length(largest_cliques(desmo_conn_graph))
length(largest_cliques(erdos_graph))
length(largest_cliques(scale_free_graph))
length(largest_cliques(graph_pa_age))
clique_num(graph_pa_age)
length(max_cliques(desmo_conn_graph, min=3))
length(max_cliques(erdos_graph, min=3))
length(max_cliques(scale_free_graph, min=3))
length(max_cliques(graph_pa_age, min=3))

pdf(file='Desmosomal_graph_weight_degree.pdf', width=10, height = 8)
{
hist_desm_degree <- hist(degree(desmo_conn_graph),freq=F,breaks=20,
                         xlab = 'degree',ylab = 'frequency', add=F,main=NA,xlim=c(0,80),ylim=c(0,0.25))
hist_desm_weight <- hist(weights_desmo,freq=F,breaks=100,
                         xlab = 'weights',ylab = 'frequency', add=F,main=NA,xlim=c(0,80),ylim=c(0,0.25))
}
dev.off()


