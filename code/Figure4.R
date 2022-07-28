#quantification of the modularity of the desmosomal connectome and random scale-free and Erdös graphs of similar statistics
#code for Figure4 of the Jasek et al desmosomal connectome paper
#Gaspar Jekely

source("code/Packages_and_Connection.R")

# Load packages
library(reticulate)
start <- Sys.time()


# read csv exported from gephi as table -----------------------------------

#this is also in the github page as desmosomal-connectome.csv
desmo_conn_table <- read.csv('data/adjacency_matrix_desmosomal_connectome_CATMAID.csv', sep = ",", header = T) 
#desmo_conn_table <- read.csv('data/desmosomal-connectome.csv', sep = ";", header = T) 

desmo_conn_matrix <- as.matrix(desmo_conn_table[2:2533],nrow=nrow(desmo_conn_table),ncol=ncol(desmo_conn_table)-1)
desmo_conn_graph <- graph_from_adjacency_matrix(desmo_conn_matrix, mode = "undirected", weighted = T,
                                                diag = TRUE, add.colnames = NULL, add.rownames = NA)

neuro_conn_table <- read.csv('data/adjacency_matrix_connectome_1042cells.csv', sep = ",", header = T) 
neuro_conn_matrix <- as.matrix(neuro_conn_table[2:1043],nrow=nrow(neuro_conn_table),ncol=ncol(neuro_conn_table)-1)
neuro_conn_graph <- graph_from_adjacency_matrix(neuro_conn_matrix, mode = "undirected", weighted = T,
                                                diag = TRUE, add.colnames = NULL, add.rownames = NA)

# graph analysis ----------------------------------------------------------

#check connected componenets
cl <- components(desmo_conn_graph)
#clustering with Leiden algorith, run with ModularityVertexPartition and resolution parameter
partition <- leiden(desmo_conn_graph, partition_type = "ModularityVertexPartition", resolution_parameter = 0.95)
partition
#calculate modularity
modularity(desmo_conn_graph,partition)

#create a new graph only from the largest connected component
desmo_conn_graph_largest = induced_subgraph(desmo_conn_graph, which(cl$membership == 1))

#sum edges between the same nodes with simplify
desmo_conn_graph_largest <- igraph::simplify(desmo_conn_graph_largest, remove.loops = FALSE,
                                             edge.attr.comb = list(weight = "sum", function(x)length(x)) )


#check connected componenets
cl_neuro <- components(neuro_conn_graph)
#create a new graph only from the largest connected component
neuro_conn_graph_largest = induced_subgraph(neuro_conn_graph, which(cl_neuro$membership == 1))

#sum edges between the same nodes with simplify
neuro_conn_graph_largest <- igraph::simplify(neuro_conn_graph_largest, remove.loops = FALSE,
                                             edge.attr.comb = list(weight = "sum", function(x)length(x)) )

partition <- leiden(desmo_conn_graph_largest, partition_type = "RBConfigurationVertexPartition",
                    resolution_parameter = 1)
#calculate modularity
system.time(modularity(desmo_conn_graph,leiden(desmo_conn_graph, resolution_parameter = 0.95)))
system.time(modularity(cluster_louvain(desmo_conn_graph)))
#use leiden method to calculate modularity

#define a modularity_leiden function
modularity_leiden <- function(graph) {
  modularity(graph,leiden(graph, resolution_parameter = 0.95))
}

# generate and analyse Erdös-Renyi graphs ---------------------------------

#sample 1000 Erdös graphs with same number of nodes and edges as the desmosomal graph 
#assign weights from the desmo graph
Erdös_graphs_1000_wg <- lapply(1:20, function(x) {
  Erdös_graph <- erdos.renyi.game(length(V(desmo_conn_graph_largest)), 
                                  length(E(desmo_conn_graph_largest)),
                                  type = "gnm",directed = FALSE,loops = FALSE)
  E(Erdös_graph)$weight <- E(desmo_conn_graph_largest)$weight
  x <- Erdös_graph
}
)

#calculate modularity scores
modularity_erdos <- lapply(Erdös_graphs_1000_wg, function(x) modularity_leiden(x))

#calculate mean distance
meandist_erdos <- lapply(Erdös_graphs_1000_wg, function(x) mean_distance(x))

#calculate diameter
diameter_erdos <- lapply(Erdös_graphs_1000_wg, function(x) diameter(x))

#degree of all nodes
degree_erdos <- lapply(Erdös_graphs_1000_wg, function(x) degree(x))

#weights of all nodes
weights_erdos <- lapply(Erdös_graphs_1000_wg, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques3_erdos <- lapply(Erdös_graphs_1000_wg, function(x) length(max_cliques(x, min=3)))

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_erdos <- lapply(Erdös_graphs_1000_wg, 
                                   function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_erdos <- lapply(Erdös_graphs_1000_wg, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))

# subsample and quantify the desmosomal graph -----------------------------


#sample subgraphs from the weighted desmosomal graph
desmo_subgraphs_1000 <- lapply(1:20, function(x)
{vids <- sample(V(desmo_conn_graph_largest),length(V(desmo_conn_graph_largest))-100)
x <- induced_subgraph(desmo_conn_graph_largest, vids, impl = "auto")}
)

#calculate modularity scores
modularity_desmo <- lapply(desmo_subgraphs_1000, function(x) modularity_leiden(x))

#calculate mean distance (same for weighted and unweighted so only for weighted)
meandist_desmo <- lapply(desmo_subgraphs_1000, function(x) mean_distance(x))

#calculate diameter
diameter_desmo <- lapply(desmo_subgraphs_1000, function(x) diameter(x))

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
V(desmo_subgraphs_1000)

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_desmo <- lapply(desmo_subgraphs_1000, 
                                   function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_desmo <- lapply(desmo_subgraphs_1000, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))

# subsample and quantify the connectome graph -----------------------------

#sample subgraphs from the weighted connectome graph
neuro_conn_subgraphs_1000 <- lapply(1:20, function(x)
{vids <- sample(V(neuro_conn_graph_largest),length(V(neuro_conn_graph_largest))-100)
x <- induced_subgraph(neuro_conn_graph_largest, vids, impl = "auto")}
)

#calculate mean distance (same for weighted and unweighted so only for weighted)
meandist_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) mean_distance(x))
mean_distance(neuro_conn_graph_largest)
diameter(neuro_conn_graph_largest)
diameter(desmo_conn_graph_largest)
#calculate diameter
diameter_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) diameter(x))
mean(as.numeric(diameter_neuro_conn))

#degree of all nodes
degree_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) degree(x))

#weights of all nodes
weights_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3 and min=4
max_cliques3_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) length(max_cliques(x, min=3)))
max(as.numeric(max_cliques3_neuro_conn))
length(max_cliques(neuro_conn_graph_largest, min=3))
length(max_cliques(neuro_conn_graph_largest, min=4))
max_cliques4_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) length(max_cliques(x, min=4)))
max(as.numeric(max_cliques4_neuro_conn))

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_neuro_conn <- lapply(neuro_conn_subgraphs_1000, 
                                        function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))

#calculate modularity score
modularity_neuro_conn <- lapply(neuro_conn_subgraphs_1000, function(x) modularity_leiden(x))

# generate and quantify scale-free graphs ---------------------------------

#sample 1000 scale-free graphs with same number of edges and nodes as the desmosomal graph and return their modularity score in a list
#assign weights from the desmo graph
sf_graphs_wg_1000 <- lapply(1:20, function(x) {
  sf_graph <- sample_fitness_pl(
    length(V(desmo_conn_graph_largest)),
    length(E(desmo_conn_graph_largest)),
    4, exponent.in = -1, loops = FALSE,
    multiple = FALSE, finite.size.correction = TRUE)
  E(sf_graph)$weight <- E(desmo_conn_graph_largest)$weight
  x=sf_graph}
)

#calculate mean distance
meandist_sf <- lapply(sf_graphs_wg_1000, function(x) mean_distance(x))
mean(as.numeric(meandist_sf))

#calculate diameter
diameter_sf <- lapply(sf_graphs_wg_1000, function(x) diameter(x))
mean(as.numeric(diameter_sf))

#degree of all nodes
degree_sf <- lapply(sf_graphs_wg_1000, function(x) degree(x))

#weights of all nodes
weights_sf <- lapply(sf_graphs_wg_1000, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques3_sf <- lapply(sf_graphs_wg_1000, function(x) length(max_cliques(x, min=3)))
max_cliques_sf4 <- lapply(sf_graphs_wg_1000, function(x) length(max_cliques(x, min=4)))

#calculate modularity score
modularity_sf <- lapply(sf_graphs_wg_1000, function(x) modularity_leiden(x))

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_sf <- lapply(sf_graphs_wg_1000, 
                                function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_sf <- lapply(sf_graphs_wg_1000, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))

# graphs with preferential attachment -------------------------------------

#sample 1000 preferential attachment graphs with same number of nodes and very similar number of edges as the desmosomal graph
#assign weights from a random sample of the desmosomal graph
pa_graphs_1000_wg <- lapply(1:20, function(x) {pa_graph <- sample_pa_age(
  length(V(desmo_conn_graph)),  pa.exp=1,aging.exp=-2,
  m = 2, aging.bin = 300,  out.dist = NULL,out.seq = NULL,
  out.pref = FALSE,directed = F,  zero.deg.appeal = 1,zero.age.appeal = 0,
  deg.coef = 1,age.coef = 1,time.window = NULL)
E(pa_graph)$weight <- sample(E(desmo_conn_graph_largest)$weight,gsize(pa_graph), replace=TRUE)
x=pa_graph
}
)

#calculate mean distance
meandist_pa <- lapply(pa_graphs_1000_wg, function(x) mean_distance(x))

#calculate diameter
diameter_pa <- lapply(pa_graphs_1000_wg, function(x) diameter(x))

#degree of all nodes
degree_pa <- lapply(pa_graphs_1000_wg, function(x) degree(x))

#weights of all nodes
weights_pa <- lapply(pa_graphs_1000_wg, function(x) strength(x))

#calculate the length of maximal cliques, run with min=3, was also run with min=4 which returns no cliques
max_cliques3_pa <- lapply(pa_graphs_1000_wg, function(x) length(max_cliques(x, min=3)))

#calculate modularity score
modularity_pa <- lapply(pa_graphs_1000_wg, function(x) modularity_leiden(x))

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_pa <- lapply(pa_graphs_1000_wg, 
                                function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_pa <- lapply(pa_graphs_1000_wg, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))


{
  #calculate node weights
  weights_desmo <- strength(desmo_conn_graph_largest, vids = V(desmo_conn_graph_largest),
                            loops = TRUE, mode = 'all')
  
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
}

# tidy the data ----------------------------------------------------------------

{
  #tidy the data for the graphs
  neuro_graphs_tb <- as_tibble(1:20) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_neuro_conn)) %>%
    mutate(diameter = unlist(diameter_neuro_conn)) %>%
    mutate(meandist = unlist(meandist_neuro_conn)) %>%
    mutate(transitivity = unlist(transitivity_centr_neuro_conn)) %>%
    mutate(three_cliques = unlist(max_cliques3_neuro_conn)) %>%
    mutate(assortativity = unlist(assort_centr_neuro_conn)) %>%
    mutate(graph_type = 'neuro')
  
  desmo_graphs_tb <- as_tibble(1:20) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_desmo)) %>%
    mutate(diameter = unlist(diameter_desmo)) %>%
    mutate(meandist = unlist(meandist_desmo)) %>%
    mutate(transitivity = unlist(transitivity_centr_desmo)) %>%
    mutate(three_cliques = unlist(max_cliques3_desmo)) %>%
    mutate(assortativity = unlist(assort_centr_desmo)) %>%
    mutate(graph_type = 'desmo')
  
  sf_graphs_tb <- as_tibble(1:20) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_sf)) %>%
    mutate(diameter = unlist(diameter_sf)) %>%
    mutate(meandist = unlist(meandist_sf)) %>%
    mutate(transitivity = unlist(transitivity_centr_sf)) %>%
    mutate(three_cliques = unlist(max_cliques3_sf)) %>%
    mutate(assortativity = unlist(assort_centr_sf)) %>%
    mutate(graph_type = 'sf')
  
  pa_graphs_tb <- as_tibble(1:20) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_pa)) %>%
    mutate(diameter = unlist(diameter_pa)) %>%
    mutate(meandist = unlist(meandist_pa)) %>%
    mutate(transitivity = unlist(transitivity_centr_pa)) %>%
    mutate(three_cliques = unlist(max_cliques3_pa)) %>%
    mutate(assortativity = unlist(assort_centr_pa)) %>%
    mutate(graph_type = 'pa')
  
  erdos_graphs_tb <- as_tibble(1:20) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_erdos)) %>%
    mutate(diameter = unlist(diameter_erdos)) %>%
    mutate(meandist = unlist(meandist_erdos)) %>%
    mutate(transitivity = unlist(transitivity_centr_erdos)) %>%
    mutate(three_cliques = unlist(max_cliques3_erdos)) %>%
    mutate(assortativity = unlist(assort_centr_erdos)) %>%
    mutate(graph_type = 'erdos')
  
  
  #join all graphs to a single tibble
  neuro_desmo_graphs_tb <- full_join(neuro_graphs_tb, desmo_graphs_tb)
  neuro_desmo_sf_graphs_tb <- full_join(neuro_desmo_graphs_tb, sf_graphs_tb)
  neuro_desmo_sf_pa_graphs_tb <- full_join(neuro_desmo_sf_graphs_tb, pa_graphs_tb)
  all_graphs_tb <- full_join(neuro_desmo_sf_pa_graphs_tb, erdos_graphs_tb)
  
  weights_desmo_tb <- as_tibble(lapply(weights_desmo, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'desmo') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "weights")
  
  weights_neuro_conn_tb <- as_tibble(lapply(weights_neuro_conn, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'neuro') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "weights")
  
  weights_sf_tb <- as_tibble(lapply(weights_sf, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'sf') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "weights")
  
  weights_pa_tb <- as_tibble(lapply(weights_pa, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'pa') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "weights")
  
  weights_erdos_tb <- as_tibble(lapply(weights_erdos, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'erdos') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "weights")
  
  weights_desmo_neuro_tb <- full_join(weights_desmo_tb, weights_neuro_conn_tb)
  weights_desmo_neuro_sf_tb <- full_join(weights_desmo_neuro_tb, weights_sf_tb)
  weights_desmo_neuro_sf_pa_tb <- full_join(weights_desmo_neuro_sf_tb, weights_pa_tb)
  weights_all_tb <- full_join(weights_desmo_neuro_sf_pa_tb, weights_erdos_tb)
  
  
  degree_desmo_tb <- as_tibble(lapply(degree_desmo, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'desmo') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "degree")
  
  degree_neuro_tb <- as_tibble(lapply(degree_neuro_conn, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'neuro') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "degree")
  
  degree_sf_tb <- as_tibble(lapply(degree_sf, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'sf') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "degree")
  
  degree_pa_tb <- as_tibble(lapply(degree_pa, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'pa') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "degree")
  
  
  degree_erdos_tb <- as_tibble(lapply(degree_erdos, function(x) x), .name_repair = "unique") %>%
    mutate(graph_type = 'erdos') %>%
    pivot_longer(cols = starts_with("..."), names_to = "replicate",
                 values_to = "degree")
  
  degree_desmo_neuro_tb <- full_join(degree_desmo_tb, degree_neuro_tb)
  degree_desmo_neuro_sf_tb <- full_join(degree_desmo_neuro_tb, degree_sf_tb)
  degree_desmo_neuro_sf_pa_tb <- full_join(degree_desmo_neuro_sf_tb, degree_pa_tb)
  degree_all_tb <- full_join(degree_desmo_neuro_sf_pa_tb, degree_erdos_tb)
}

# plot degree distributions ------------------------------------

{
  
  weights_all_tb %>%
    ggplot() +
    geom_freqpoly(aes(x = weights, after_stat(density), 
                      color = graph_type, alpha = replicate), 
                  binwidth = 0.2) +
    scale_x_log10() +
    theme_minimal() +
    xlab('degree') +
    theme()
  
  
  degree_all_tb %>%
    ggplot() +
    geom_freqpoly(aes(x = degree, after_stat(density), 
                      color = graph_type, alpha = replicate), 
                  binwidth = 0.2) +
    scale_x_log10() +
    theme_minimal() +
    xlab('degree') +
    theme()  
  
  ggplot() +
    geom_freqpoly(data = as_tibble(unlist(degree_desmo)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#0072B2", 
                  alpha = 0.6,  size = 1, linetype = 2) +
    geom_freqpoly(data = as_tibble(unlist(degree_neuro_conn)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#E69F00", 
                  alpha = 0.7, size = 1, linetype = 1) +
    geom_freqpoly(data = as_tibble(unlist(degree_sf)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#000000", 
                  alpha = 0.6, size = 1, linetype = 1) +
    geom_freqpoly(data = as_tibble(unlist(degree_pa)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "grey90", 
                  alpha = 1, size = 2, linetype = 4) +
    geom_freqpoly(data = as_tibble(unlist(degree_erdos)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, color = "#D55E00", 
                  alpha = 0.8, size = 1, linetype = 5)  +
    scale_x_log10() +
    theme_minimal() +
    xlab('degree') +
    theme(legend.position = c(0.2, .8))
  ggsave('pictures/degree.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot weights --------------------------------------------------------

weights_all_tb
{
  ggplot() +
    geom_freqpoly(data = as_tibble(unlist(weights_desmo)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#0072B2", 
                  alpha = 0.6,  size = 1, linetype = 2) +
    geom_freqpoly(data = as_tibble(unlist(weights_neuro_conn)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#E69F00", 
                  alpha = 0.7, size = 1, linetype = 1) +
    geom_freqpoly(data = as_tibble(unlist(weights_sf)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "#000000", 
                  alpha = 0.6, size = 1, linetype = 1) +
    geom_freqpoly(data = as_tibble(unlist(weights_pa)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, 
                  color = "grey90", 
                  alpha = 1, size = 2, linetype = 4) +
    geom_freqpoly(data = as_tibble(unlist(weights_erdos)),
                  aes(x = value, after_stat(density)), binwidth = 0.2, color = "#D55E00", 
                  alpha = 0.8, size = 1, linetype = 5)  +
    scale_x_log10() +
    theme_minimal() +
    xlab('weight')
  ggsave('pictures/weights.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot modularity ------------------------------------

{
  modul_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = modularity, fill = graph_type), 
                   binwidth = 0.01, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type"))
  
  ggsave('pictures/modularity.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot histograms of max clique scores ------------------------------------

{
  max3_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = three_cliques, fill = graph_type), 
                   binwidth = 0.1, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('# of 3-cliques')
  
  ggsave('pictures/3cliques.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot transitivity --------------------------------------------------------

{
  trans_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = transitivity, fill = graph_type), 
                   binwidth = 0.1, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('transitivity')
  
  ggsave('pictures/transitivity.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot assortativity --------------------------------------------------------

{
  assort_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = assortativity, fill = graph_type), 
                   binwidth = 0.1, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('assortativity')
  
  
  ggsave('pictures/assortativity.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot diameter --------------------------------------------------------

{
  diameter_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = diameter, fill = graph_type), 
                   binwidth = 0.01, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('diameter')
  
  ggsave('pictures/diameter.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot meandist --------------------------------------------------------

{
  meandist_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = meandist, fill = graph_type), 
                   binwidth = 0.01, 
                   alpha = 0.6, color = "grey20", size = 0.1) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('meandist')
  
  ggsave('pictures/meandist.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# create multi-panel figure -----------------------------------------------
{
  Fig4_degree <- ggdraw() + draw_image(readPNG("pictures/degree.png"))
  
  Fig4_weights <- ggdraw() + draw_image(readPNG("pictures/weights.png"))
  
  Fig4_trans <- ggdraw() + draw_image(readPNG("pictures/transitivity.png")) +
    draw_label("desmo", x = 0.72, y = 0.71, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.91, y = 0.63, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.5, y = 0.75, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.66, y = 0.5, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.3, y = 0.49, fontfamily = "sans", size = 8)
  
  Fig4_diameter <- ggdraw() + draw_image(readPNG("pictures/diameter.png")) +
    draw_label("desmo", x = 0.88, y = 0.63, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.7, y = 0.55, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.24, y = 0.84, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.37, y = 0.58, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.31, y = 0.72, fontfamily = "sans", size = 8)
  
  Fig4_meandist <- ggdraw() + draw_image(readPNG("pictures/meandist.png")) +
    draw_label("desmo", x = 0.88, y = 0.64, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.25, y = 0.78, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.545, y = 0.64, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.72, y = 0.48, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.6, y = 0.74, fontfamily = "sans", size = 8)
  
  Fig4_3cl <- ggdraw() + draw_image(readPNG("pictures/3cliques.png")) +
    draw_label("desmo", x = 0.76, y = 0.61, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.91, y = 0.8, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.48, y = 0.55, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.62, y = 0.5, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.29, y = 0.75, fontfamily = "sans", size = 8)
  
  Fig4_modul <- ggdraw() + draw_image(readPNG("pictures/modularity.png")) +
    draw_label("desmo", x = 0.9, y = 0.67, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.49, y = 0.6, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.26, y = 0.5, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.76, y = 0.5, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.38, y = 0.75, fontfamily = "sans", size = 8)
  
  Fig4_assort <- ggdraw() + draw_image(readPNG("pictures/assortativity.png")) +
    draw_label("desmo", x = 0.24, y = 0.84, fontfamily = "sans", size = 8) +
    draw_label("neuro", x = 0.6, y = 0.55, fontfamily = "sans", size = 8) +
    draw_label("sf", x = 0.46, y = 0.64, fontfamily = "sans", size = 8) +
    draw_label("pa", x = 0.88, y = 0.63, fontfamily = "sans", size = 8) +
    draw_label("Erdös", x = 0.41, y = 0.72, fontfamily = "sans", size = 8)
  
  layout <- "
ABCD
EFGH"
  
  Fig4 <- Fig4_degree + Fig4_weights + modul_plot + diameter_plot + 
    meandist_plot + trans_plot + max3_plot + assort_plot +
    plot_layout(design = layout, heights = c(1, 1),
                guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 12, face='plain'))
  
  Fig4 <- plot_grid(Fig4_degree, Fig4_weights, Fig4_modul, Fig4_diameter,
                    Fig4_meandist, Fig4_trans, Fig4_3cl, Fig4_assort,
                    ncol=4,
                    rel_heights = c(1, 1,1,1),
                    labels=c("A", "B","C","D","E","F","G","H"),
                    label_size = 12, label_y = 1, label_x = -0.04,
                    label_fontfamily = "sans", label_fontface = "plain") + 
    theme(plot.margin = unit(c(1, 1, 1, 3), units = "pt"))
  
  ggsave("figures/Figure_4.pdf", limitsize = FALSE, 
         units = c("px"), Fig4, width = 2600, height = 1400)
}
ggsave("figures/Figure_4.png", limitsize = FALSE, 
       units = c("px"), Fig4, width = 2600, height = 1400)



