#quantification of the modularity of the desmosomal connectome and random scale-free and Erdos graphs of similar statistics
#code for Figure4 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely 2021 March

#set working directory

setwd('/workdir/')

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

# Load packages
library(igraph)
#https://rdrr.io/cran/igraph/man/
library(leiden)
#https://github.com/TomKellyGenetics/leiden

#check connected componenets
cl <- components(desmo_conn_graph)
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
sum(E(c)$weight)

#convert into binary matrix
desmo_conn_matrix_bi<-as.matrix((desmo_conn_matrix>0)+0)
desmo_conn_graph_bi <- graph_from_adjacency_matrix(desmo_conn_matrix_bi, mode = "undirected", weighted = T,
                                                   diag = TRUE, add.colnames = NULL, add.rownames = NA)

#create a new graph only from the largest connected component
desmo_conn_graph_bi_largest = induced_subgraph(desmo_conn_graph_bi, which(cl$membership == 1))

partition <- leiden(desmo_conn_graph_largest, partition_type = "RBConfigurationVertexPartition",
                    resolution_parameter = 1)
max(partition)
#calculate modularity


system.time(modularity(desmo_conn_graph,leiden(desmo_conn_graph, resolution_parameter = 0.95)))
system.time(modularity(cluster_louvain(desmo_conn_graph)))
#will use louvain clustering to calculate modularity - for revision recalculated with Leiden

#define a modularity_leiden function
modularity_leiden <- function(graph) {
  modularity(graph,leiden(graph, resolution_parameter = 0.95))
}

#######################################
#generate and analyse Erdos-Renyi graphs

#sample 1000 Erdos graphs with same number of nodes and edges as the desmosomal graph 
erdos_graphs_1000 <- lapply(1:10, function(x) 
  x=erdos.renyi.game(length(V(desmo_conn_graph_largest)), 
                     length(E(desmo_conn_graph_largest)),
                     type = "gnm",directed = FALSE,loops = FALSE))

#do the same but assing weights from the desmo graph
erdos_graphs_1000_wg <- lapply(1:10, function(x) {
  erdos_graph <- erdos.renyi.game(length(V(desmo_conn_graph_largest)), 
                                  length(E(desmo_conn_graph_largest)),
                                  type = "gnm",directed = FALSE,loops = FALSE)
  E(erdos_graph)$weight <- E(desmo_conn_graph_largest)$weight
  x <- erdos_graph
}
)

#calculate modularity scores
modularity_erdos <- lapply(erdos_graphs_1000, function(x) modularity_leiden(x))
modularity_erdos_wg <- lapply(erdos_graphs_1000_wg, function(x) modularity_leiden(x))

#calculate mean distance
meandist_erdos <- lapply(erdos_graphs_1000, function(x) mean_distance(x))
mean(as.numeric(meandist_erdos))

#calculate diameter
diameter_erdos <- lapply(erdos_graphs_1000, function(x) diameter(x))
mean(as.numeric(diameter_erdos))

#degree of all nodes
degree_erdos <- lapply(erdos_graphs_1000, function(x) degree(x))

#weights of all nodes
weights_erdos <- lapply(erdos_graphs_1000_wg, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques_erdos <- lapply(erdos_graphs_1000, function(x) length(max_cliques(x, min=3)))
max(as.numeric(max_cliques_erdos))

#calculate the global transitivity (clustering coefficient) 
transitivity_erdos <- lapply(erdos_graphs_1000, function(x) transitivity(x))

#calculate eigen vector centrality 
eigen_centr_erdos <- lapply(erdos_graphs_1000, function(x) eigen_centrality(x, directed = FALSE, weights = NA, scale = TRUE))


################################
#subsample and quantify the desmosomal graph

#sample subgraphs from the weighted desmosomal graph
desmo_subgraphs_1000 <- lapply(1:10, function(x)
{vids <- sample(V(desmo_conn_graph_largest),length(V(desmo_conn_graph_largest))-100)
x <- induced_subgraph(desmo_conn_graph_largest, vids, impl = "auto")}
)

#sample subgraphs from the binarised desmosomal graph
desmo_bi_subgraphs_1000 <- lapply(1:10, function(x)
{vids <- sample(V(desmo_conn_graph_bi_largest),length(V(desmo_conn_graph_largest))-100)
x <- induced_subgraph(desmo_conn_graph_bi_largest, vids, impl = "auto")}
)

#calculate modularity scores
modularity_desmo <- lapply(desmo_subgraphs_1000, function(x) modularity_leiden(x))
modularity_desmo_bi <- lapply(desmo_bi_subgraphs_1000, function(x) modularity_leiden(x))

#calculate mean distance (same for weighted and unweighted so only for weighted)
meandist_desmo <- lapply(desmo_subgraphs_1000, function(x) mean_distance(x))
mean_distance(desmo_conn_graph_largest)

#calculate diameter
diameter_desmo <- lapply(desmo_subgraphs_1000, function(x) diameter(x))
mean(as.numeric(diameter_desmo))

#degree of all nodes
degree_desmo <- lapply(desmo_subgraphs_1000, function(x) degree(x))

#weights of all nodes
weights_desmo <- lapply(desmo_subgraphs_1000, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3 and min=4
max_cliques3_desmo <- lapply(desmo_subgraphs_1000, function(x) length(max_cliques(x, min=3)))
max(as.numeric(max_cliques3_desmo))
length(max_cliques(desmo_conn_graph_largest, min=3))
length(max_cliques(desmo_conn_graph_largest, min=4))
max_cliques4_desmo <- lapply(desmo_subgraphs_1000, function(x) length(max_cliques(x, min=4)))
max(as.numeric(max_cliques4_desmo))

#calculate the global transitivity (clustering coefficient) 
transitivity_desmo <- lapply(desmo_subgraphs_1000, function(x) transitivity(x))
transitivity_desmo_bi <- lapply(desmo_bi_subgraphs_1000, function(x) transitivity(x))

#calculate eigen vector centrality 
eigen_centr_desmo <- lapply(desmo_subgraphs_1000, function(x) eigen_centrality(x, directed = FALSE, weights = NA, scale = TRUE))



#######################################
#generate and quantify scale-free graphs

#sample 1000 scale-free graphs with same number of edges and nodes as the desmosomal graph and return their modularity score in a list
sf_graphs_1000 <- lapply(1:10, function(x) 
  x=sample_fitness_pl(
    length(V(desmo_conn_graph_largest)),
    length(E(desmo_conn_graph_largest)),
    4,
    exponent.in = -1,
    loops = FALSE,
    multiple = FALSE,
    finite.size.correction = TRUE
  ))

#do the same but assing weights from the desmo graph
sf_graphs_wg_1000 <- lapply(1:10, function(x) {
  sf_graph <- sample_fitness_pl(
    length(V(desmo_conn_graph_largest)),
    length(E(desmo_conn_graph_largest)),
    4, exponent.in = -1, loops = FALSE,
    multiple = FALSE, finite.size.correction = TRUE)
  E(sf_graph)$weight <- E(desmo_conn_graph_largest)$weight
  x=sf_graph}
)

#calculate mean distance
meandist_sf <- lapply(sf_graphs_1000, function(x) mean_distance(x))
mean(as.numeric(meandist_sf))

#calculate diameter
diameter_sf <- lapply(sf_graphs_1000, function(x) diameter(x))
mean(as.numeric(diameter_sf))

#degree of all nodes
degree_sf <- lapply(sf_graphs_1000, function(x) degree(x))

#weights of all nodes
weights_sf <- lapply(sf_graphs_wg_1000, function(x) strength(x))


#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques_sf <- lapply(sf_graphs_1000, function(x) length(max_cliques(x, min=3)))
max(as.numeric(max_cliques_sf))
max_cliques_sf4 <- lapply(sf_graphs_1000, function(x) length(max_cliques(x, min=4)))
max(as.numeric(max_cliques_sf))

#calculate the global transitivity (clustering coefficient) 
transitivity_sf <- lapply(sf_graphs_1000, function(x) transitivity(x))
mean(as.numeric(transitivity_sf))

#calculate modularity score
modularity_sf <- lapply(sf_graphs_1000, function(x) modularity_leiden(x))
modularity_sf_wg <- lapply(sf_graphs_wg_1000, function(x) modularity_leiden(x))

#calculate eigen vector centrality 
eigen_centr_sf <- lapply(sf_graphs_1000, function(x) eigen_centrality(x, directed = FALSE, weights = NA, scale = TRUE))


#####################################
#graphs with preferential attachment

#sample 1000 preferential attachment graphs with same number of nodes and very similar number of edges as the desmosomal graph
pa_graphs_1000 <- lapply(1:10, function(x) x=(sample_pa_age(
  length(V(desmo_conn_graph)),  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL)
))

#do the same and assign weights from a random sample of the desmosomal graph
pa_graphs_1000_wg <- lapply(1:10, function(x) {pa_graph <- sample_pa_age(
  length(V(desmo_conn_graph)),  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL)
E(pa_graph)$weight <- sample(E(desmo_conn_graph_largest)$weight,gsize(pa_graph), replace=TRUE)
x=pa_graph
}
)

#calculate mean distance
meandist_pa <- lapply(pa_graphs_1000, function(x) mean_distance(x))
mean(as.numeric(meandist_pa))

#calculate diameter
diameter_pa <- lapply(pa_graphs_1000, function(x) diameter(x))
mean(as.numeric(diameter_pa))

#degree of all nodes
degree_pa <- lapply(pa_graphs_1000, function(x) degree(x))

#weights of all nodes
weights_pa <- lapply(pa_graphs_1000_wg, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques_pa <- lapply(pa_graphs_1000, function(x) length(max_cliques(x, min=3)))
max(as.numeric(max_cliques_pa))

#calculate the global transitivity (clustering coefficient) 
transitivity_pa <- lapply(pa_graphs_1000, function(x) transitivity(x))

#calculate modularity score
modularity_pa <- lapply(pa_graphs_1000, function(x) modularity_leiden(x))
modularity_pa_wg <- lapply(pa_graphs_1000_wg, function(x) modularity_leiden(x))
mean(as.numeric(modularity_pa))

#calculate eigen vector centrality 
eigen_centr_pa <- lapply(pa_graphs_1000, function(x) eigen_centrality(x, directed = FALSE, weights = NA, scale = TRUE))
eigen_centr_pa_wg <- lapply(pa_graphs_1000_wg, function(x) eigen_centrality(x, directed = FALSE, weights = NA, scale = TRUE))


#############################################
#plot histograms of modularity scores


{
  hist_desmo_subgraphs <- hist(as.numeric(modularity_desmo),
                               xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_desmo_bi_subgraphs <- hist(as.numeric(modularity_desmo_bi),
                                  xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_erdos <- hist(as.numeric(modularity_erdos),
                     xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_erdos_weighted <- hist(as.numeric(modularity_erdos_wg),
                              xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_sf <- hist(as.numeric(modularity_sf),
                  xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_sf_weighted <- hist(as.numeric(modularity_sf_wg),
                           xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_pa <- hist(as.numeric(modularity_pa),
                  xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
  hist_pa_weighted <- hist(as.numeric(modularity_pa_wg),
                           xlim=c(0.5,0.9),ylim=c(0,400),ylab = 'count',xlab='modularity', main=NA)
}

pdf(file='Desmosomal_etc_graph_wg_modularity.pdf', width=8, height = 8)
{
  par(mar=c(6,6,2,2)) 
  plot(hist_desmo_subgraphs,add=F,xlim=c(0.7,0.9),ylim=c(0,250),ylab = '',xlab='', cex.axis=2,main=NA)
  plot(hist_erdos_weighted,add=T)
  plot(hist_sf_weighted,add=T, col = hcl.colors(1,palette = 'Blues',alpha=0.4), border = hcl.colors(1,palette = 'Blues',alpha=0.4))
  plot(hist_pa_weighted,add=T, col = hcl.colors(1,palette = 'Reds',alpha=0.4), border = hcl.colors(1,palette = 'Reds',alpha=0.5))
  title(xlab='modularity, weighted', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()

pdf(file='Desmosomal_etc_graph_unwg_modularity.pdf', width=8, height = 8)
{
  par(mar=c(6,6,2,2)) 
  plot(hist_desmo_bi_subgraphs,add=F,xlim=c(0.5,0.9),ylim=c(0,400),ylab = '',xlab='', cex.axis=2,main=NA)
  plot(hist_erdos,add=T)
  plot(hist_sf,add=T, col = hcl.colors(1,palette = 'Blues',alpha=0.4), border = hcl.colors(1,palette = 'Blues',alpha=0.4))
  plot(hist_pa,add=T,col = hcl.colors(1,palette = 'Reds',alpha=0.5), border = hcl.colors(1,palette = 'Reds',alpha=0.5))
  title(xlab='modularity, unweighted', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()


#calculate median modularity values
median(as.numeric(modularity_erdos))
median(as.numeric(modularity_erdos_wg))
median(as.numeric(modularity_desmo))
median(as.numeric(modularity_desmo_bi))
median(as.numeric(modularity_sf))
median(as.numeric(modularity_sf_wg))
median(as.numeric(modularity_pa))
median(as.numeric(modularity_pa_wg))
modularity_leiden(desmo_conn_graph)

#check package versions
packageVersion("igraph") 
packageVersion("leiden") 


#############################################
#plot histograms of transitivity scores

{
  hist_desmo_transitivity <- hist(as.numeric(transitivity_desmo),
                                  xlim=c(0,0.09),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA)
  hist_erdos_transitivity <- hist(as.numeric(transitivity_erdos),
                                  xlim=c(0,0.09),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA,col='red')
  hist_sf_transitivity <- hist(as.numeric(transitivity_sf),
                               xlim=c(0,0.09),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA, col='cyan')
  hist_pa_transitivity <- hist(as.numeric(transitivity_pa),
                               xlim=c(0,0.09),ylim=c(0,400),ylab = 'count',xlab='global transitivity', main=NA, col='cyan')
}

pdf(file='Desmosomal_Erdos_and_sf_graph_transitivity.pdf', width=8, height = 8)
{
  par(mar=c(6,6,2,2)) 
  plot(hist_desmo_transitivity,add=F,xlim=c(0,0.09),ylim=c(0,400),ylab = '',xlab='', cex.axis=2, main=NA)
  plot(hist_erdos_transitivity,add=T)
  plot(hist_sf_transitivity,add=T, col = hcl.colors(1,palette = 'Blues 2',alpha=0.4), border = hcl.colors(1,palette = 'Blues 2',alpha=1))
  plot(hist_pa_transitivity,add=T, , col = hcl.colors(1,palette = 'Reds',alpha=0.4), border = hcl.colors(1,palette = 'Reds',alpha=0.4))
  title(xlab='global transitivity', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()

median(as.numeric(transitivity_desmo))
median(as.numeric(transitivity_erdos))
median(as.numeric(transitivity_sf))
median(as.numeric(transitivity_pa))


################################################
#plot histograms of max clique scores

{
  hist_desmo_cliques3 <- hist(as.numeric(max_cliques3_desmo),
                              xlim=c(0,600),ylim=c(0,300),ylab = 'count',xlab='cliques', main=NA)
  hist_desmo_cliques4 <- hist(as.numeric(max_cliques4_desmo),
                              xlim=c(0,600),ylim=c(0,300),ylab = 'count',xlab='cliques', main=NA)
  hist_erdos_cliques3 <- hist(as.numeric(max_cliques_erdos),
                              xlim=c(0,600),ylim=c(0,300),ylab = 'count',xlab='cliques', main=NA)
  hist_sf_cliques3 <- hist(as.numeric(max_cliques_sf),
                           xlim=c(0,600),ylim=c(0,300),ylab = 'count',xlab='cliques', main=NA)
  hist_pa_cliques3 <- hist(as.numeric(max_cliques_pa),
                           xlim=c(0,600),ylim=c(0,300),ylab = 'count',xlab='cliques', main=NA)
}

pdf(file='Desmosomal_Erdos_sf_pa_graphs_3cliques.pdf', width=8, height = 8)
{
  par(mar=c(6,6,2,2)) 
  plot(hist_desmo_cliques3,add=F,xlim=c(0,600),ylim=c(0,400),ylab = '',xlab='',
       cex.axis=2, main=NA)
  plot(hist_erdos_cliques3,add=T, border = hcl.colors(1,palette = 'Grays',alpha=0.3))
  plot(hist_sf_cliques3,add=T,col = hcl.colors(1,palette = 'Blues 2',alpha=1), border = hcl.colors(1,palette = 'Blues 2',alpha=1))
  #plot(hist_desmo_cliques4,add=T,border='red')
  plot(hist_pa_cliques3,add=T, col = hcl.colors(1,palette = 'Reds',alpha=0.4), border = hcl.colors(1,palette = 'Reds',alpha=1))
  title(xlab='# of 3-cliques', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()


#calculate node weights
weights_desmo <- strength(desmo_conn_graph_largest, vids = V(desmo_conn_graph_largest),
                          loops = TRUE)

#other network parameters
mean(weights_desmo)
diameter(desmo_conn_graph_largest)
mean_distance(desmo_conn_graph)
length(which((coreness(desmo_conn_graph,mode='all'))==5))
max(weights_desmo)
sort(degree(desmo_conn_graph), decreasing=T)
count_motifs(desmo_conn_graph, 4)
max(count_triangles(desmo_conn_graph))
(triangles(desmo_conn_graph))

#triangles
max(count_triangles(desmo_conn_graph))
count_triangles(desmo_conn_graph, vids=V(desmo_conn_graph)[223])
V(desmo_conn_graph)[223]

#cliques
length(largest_cliques(desmo_conn_graph))
clique_num(desmo_conn_graph)
length(max_cliques(desmo_conn_graph, min=3))

#########################################
#plot degree and weight distributions


pdf(file='Desmosomal_graph_weight_degree.pdf', width=8, height = 8)
{
  par(mar=c(6,6,2,2)) 
  hist_desm_degree <- hist(degree(desmo_conn_graph_largest),freq=F,breaks=20,
                           xlab = '',ylab = '', cex.axis=2, add=F,main=NA,xlim=c(0,80),
                           ylim=c(0,0.45))
  title(xlab='degree', ylab = 'frequency', mgp = c(4,1, 0),cex.lab=3, font.lab=2)
  par(mar=c(6,6,2,2)) 
  hist_desm_weight <- hist(weights_desmo,freq=F,breaks=60,
                           xlab = '',ylab = '', cex.axis=2, add=F,main=NA,xlim=c(0,80),
                           ylim=c(0,0.45))
  title(xlab='weights', ylab = 'frequency', mgp = c(4,1, 0),cex.lab=3, font.lab=2)
}
dev.off()


pdf(file='Degree_distr.pdf', width=8, height = 8)
{ hist_degree_desmo <- hist(unlist(degree_desmo), freq=FALSE,
                            breaks=c(0:100))
  hist_degree_sf <- hist(unlist(degree_sf), freq=FALSE,
                         breaks=c(0:150))
  hist_degree_pa <- hist(unlist(degree_pa), freq=FALSE,
                         breaks=c(0:150))
  hist_degree_erdos <- hist(unlist(degree_erdos), freq=FALSE,
                            breaks=c(0:150))
  hist_degree_desmo
  par(mar=c(6,6,2,2)) 
  plot(hist_degree_desmo$density,xlim=c(0,40),ylim=c(0,0.35),
       main=NA, cex.axis=2, xlab='', ylab='',col='grey30',lwd=5)
  lines(hist_degree_sf$density, 
        col = 'blue',lty=1,lwd=8)
  lines(hist_degree_pa$density,
        col = 'red', lwd=5)
  lines(hist_degree_erdos$density, 
        col = hcl.colors(1,palette = 'Grays',alpha=0.5) ,lty=6,lwd=7)
  title(xlab='degree', ylab = 'frequency', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
  legend("topright", inset=.02, title="",
         c("desmo", "sf", "pa", "Erdös"), lwd=c(12,5,6,7), cex=2, lty=c(3,1,1,6),
         col=c('black','blue','red','grey40'),bty='n',
         x.intersp=1, ncol=1)
}
dev.off()


pdf(file='Weight_distr.pdf', width=8, height = 8)
{  hist_weight_desmo <- hist(unlist(weights_desmo), freq=FALSE,
                             breaks=c(0:100))
  hist_weight_sf <- hist(unlist(weights_sf), freq=FALSE,
                         breaks=c(0:150))
  hist_weight_pa <- hist(unlist(weights_pa), freq=FALSE,
                         breaks=c(0:150))
  hist_weight_erdos <- hist(unlist(weights_erdos), freq=FALSE,
                            breaks=c(0:150))
  par(mar=c(6,6,2,2)) 
  plot(hist_weight_desmo$density,xlim=c(0,40),ylim=c(0,0.13),
       main=NA, cex.axis=2, xlab='', ylab='',col='grey30',lwd=5)
  lines(hist_weight_sf$density, 
        col = 'blue',lty=1,lwd=8)
  lines(hist_weight_pa$density,
        col = 'red', lwd=5)
  lines(hist_weight_erdos$density, 
        col = hcl.colors(1,palette = 'Grays',alpha=0.5) ,lty=6,lwd=7)
  title(xlab='weight', ylab = 'frequency', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
  legend("topright", inset=.02, title="",
         c("desmo","sf", "pa", "Erdös"), lwd=c(12,5,6,7), cex=2, lty=c(3,1,1,6),
         col=c('black','blue','red','grey40'),bty='n',
         x.intersp=1, ncol=1)
}
dev.off()


#####################################################
#eigenvalues

pdf(file='Eigenvalues.pdf', width=8, height = 8)
{
  hist_eigen_desmo <- hist(as.numeric(lapply(eigen_centr_desmo, function(x) x[]$value)))
  hist_eigen_sf <- hist(as.numeric(lapply(eigen_centr_sf, function(x) x[]$value)))
  hist_eigen_pa <- hist(as.numeric(lapply(eigen_centr_pa, function(x) x[]$value)))
  hist_eigen_erdos <- hist(as.numeric(lapply(eigen_centr_erdos, function(x) x[]$value)))
  par(mar=c(6,6,2,2)) 
  plot(hist_eigen_desmo,add=F,xlim=c(4,10),ylim=c(0,300),
       main=NA, cex.axis=2, xlab='', ylab='')
  plot(hist_eigen_sf,add=T, col = hcl.colors(1,palette = 'Blues',alpha=0.4))
  plot(hist_eigen_pa,add=T,col = hcl.colors(1,palette = 'Reds',alpha=0.5))
  plot(hist_eigen_erdos,add=T, border = hcl.colors(1,palette = 'Grays',alpha=0.5) )
  title(xlab='Eigenvalue', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()

#####################################################
#eigen centrality values

pdf(file='Meandistance_values.pdf', width=8, height = 8)
{
  hist_meandist_desmo <- hist(as.numeric(meandist_desmo))
  hist_meandist_sf <- hist(as.numeric(meandist_sf))
  hist_meandist_pa <- hist(as.numeric(meandist_pa))
  hist_meandist_erdos <- hist(as.numeric(meandist_erdos))
  par(mar=c(6,6,2,2)) 
  plot(hist_meandist_desmo,add=F,xlim=c(5,10),ylim=c(0,300),
       main=NA, cex.axis=2, xlab='', ylab='')
  plot(hist_meandist_sf,add=T, col = hcl.colors(1,palette = 'Blues',alpha=0.4),border = hcl.colors(1,palette = 'Blues',alpha=0.4))
  plot(hist_meandist_pa,add=T,col = hcl.colors(1,palette = 'Reds',alpha=0.5))
  plot(hist_meandist_erdos,add=T, border = hcl.colors(1,palette = 'Grays',alpha=0.5) )
  title(xlab='mean distance', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()


#####################################################
#plot diameter

pdf(file='diameter_values.pdf', width=8, height = 8)
{
  hist_diameter_desmo <- hist(as.numeric(diameter_desmo))
  hist_diameter_sf <- hist(as.numeric(diameter_sf))
  hist_diameter_pa <- hist(as.numeric(diameter_pa))
  hist_diameter_erdos <- hist(as.numeric(diameter_erdos))
  par(mar=c(6,6,2,2)) 
  plot(hist_diameter_desmo,add=F,xlim=c(0,90),ylim=c(0,500),
       main=NA, cex.axis=2, xlab='', ylab='')
  plot(hist_diameter_sf,add=T, col = hcl.colors(1,palette = 'Blues',alpha=0.4),border = hcl.colors(1,palette = 'Blues',alpha=0.4))
  plot(hist_diameter_pa,add=T,col = hcl.colors(1,palette = 'Reds',alpha=0.5),border  = hcl.colors(1,palette = 'Reds',alpha=0.5))
  plot(hist_diameter_erdos,add=T, border = hcl.colors(1,palette = 'Grays',alpha=0.5) )
  title(xlab='diameter', ylab = 'count', mgp = c(4,1, 0),cex.lab=3, font.lab=2) 
}
dev.off()


