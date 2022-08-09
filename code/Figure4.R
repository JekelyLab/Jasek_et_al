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
Erdös_graphs_1000_wg <- lapply(1:1000, function(x) {
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
desmo_subgraphs_1000 <- lapply(1:1000, function(x)
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

#calculate transitivity centrality (clustering coefficient) 
transitivity_centr_desmo <- lapply(desmo_subgraphs_1000, 
                                   function(x) transitivity(x, type = "average", weights = NULL))

#calculate assortativity 
assort_centr_desmo <- lapply(desmo_subgraphs_1000, function(x) 
  assortativity(x, directed = FALSE, types1 = V(x), types2 = NULL))

# subsample and quantify the connectome graph -----------------------------

{
#sample subgraphs from the weighted connectome graph
neuro_conn_subgraphs_1000 <- lapply(1:1000, function(x)
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
}

# generate and quantify scale-free graphs ---------------------------------
{
#sample 1000 scale-free graphs with same number of edges and nodes as the desmosomal graph and return their modularity score in a list
#assign weights from the desmo graph
sf_graphs_wg_1000 <- lapply(1:1000, function(x) {
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
}
# graphs with preferential attachment -------------------------------------
{
#sample 1000 preferential attachment graphs with same number of nodes and very similar number of edges as the desmosomal graph
#assign weights from a random sample of the desmosomal graph
pa_graphs_1000_wg <- lapply(1:1000, function(x) {pa_graph <- sample_pa_age(
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
}

beep()
# tidy the data ----------------------------------------------------------------

{
#tidy the data for the graphs
neuro_graphs_tb <- as_tibble(1:1000) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_neuro_conn)) %>%
    mutate(diameter = unlist(diameter_neuro_conn)) %>%
    mutate(meandist = unlist(meandist_neuro_conn)) %>%
    mutate(transitivity = unlist(transitivity_centr_neuro_conn)) %>%
    mutate(three_cliques = unlist(max_cliques3_neuro_conn)) %>%
    mutate(assortativity = unlist(assort_centr_neuro_conn)) %>%
    mutate(graph_type = 'neuro')
  
desmo_graphs_tb <- as_tibble(1:1000) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_desmo)) %>%
    mutate(diameter = unlist(diameter_desmo)) %>%
    mutate(meandist = unlist(meandist_desmo)) %>%
    mutate(transitivity = unlist(transitivity_centr_desmo)) %>%
    mutate(three_cliques = unlist(max_cliques3_desmo)) %>%
    mutate(assortativity = unlist(assort_centr_desmo)) %>%
    mutate(graph_type = 'desmo')
  
sf_graphs_tb <- as_tibble(1:1000) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_sf)) %>%
    mutate(diameter = unlist(diameter_sf)) %>%
    mutate(meandist = unlist(meandist_sf)) %>%
    mutate(transitivity = unlist(transitivity_centr_sf)) %>%
    mutate(three_cliques = unlist(max_cliques3_sf)) %>%
    mutate(assortativity = unlist(assort_centr_sf)) %>%
    mutate(graph_type = 'sf')
  
pa_graphs_tb <- as_tibble(1:1000) %>%
    rename('graph_number'=value) %>%
    mutate(modularity = unlist(modularity_pa)) %>%
    mutate(diameter = unlist(diameter_pa)) %>%
    mutate(meandist = unlist(meandist_pa)) %>%
    mutate(transitivity = unlist(transitivity_centr_pa)) %>%
    mutate(three_cliques = unlist(max_cliques3_pa)) %>%
    mutate(assortativity = unlist(assort_centr_pa)) %>%
    mutate(graph_type = 'pa')
  
erdos_graphs_tb <- as_tibble(1:1000) %>%
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

# plot edge weight and node degree distributions ------------------------------------

{
  degree_plot <- degree_all_tb %>%
    ggplot() +
    geom_freqpoly(aes(x = degree, after_stat(density), 
                      color = graph_type, alpha = replicate), 
                  binwidth = 0.21) +
    scale_x_sqrt(breaks = c(1, 10, 50, 100)) +
    theme_minimal() +
    labs(x = 'node degree', color = 'graph type') +
    theme_minimal() +
    scale_color_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                  'grey70', '#000000'),
                       breaks = c('desmo', 'erdos', 'neuro', 
                                  'pa', 'sf')) +
    scale_linetype_manual(values = c(1, 2, 3,
                                     4, 5),
                          breaks = c('desmo', 'erdos', 'neuro', 
                                     'pa', 'sf')) +
    scale_alpha_manual(values = c(250:1250)/100000) +
    guides(color = "none", alpha = "none")
  
ggsave('pictures/degree.png', degree_plot, limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
  
weights_plot <- weights_all_tb %>%
    ggplot() +
    geom_freqpoly(aes(x = weights, after_stat(density), 
                      color = graph_type, alpha = replicate), 
                  binwidth = 0.21) +
    scale_x_sqrt(breaks = c(1, 10, 100, 400)) +
    theme_minimal() +
    labs(x = 'edge weight', color = 'graph type') +
    theme_minimal() +
    scale_color_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey70', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    scale_linetype_manual(values = c(1, 2, 3,
                                     4, 5),
                          breaks = c('desmo', 'erdos', 'neuro', 
                                     'pa', 'sf')) +
    scale_alpha_manual(values = c(250:1250)/100000) +
    guides(alpha = "none")

ggsave('pictures/weights.png', weights_plot, limitsize = FALSE, 
       units = c("px"), width = 1100, height = 800, bg = 'white')

}

# plot modularity ------------------------------------

{
    
modul_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = modularity, fill = graph_type), 
                   binwidth = 0.005, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    scale_x_continuous() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = "none")

modul_plot  
ggsave('pictures/modularity.png', modul_plot, limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot histograms of max clique scores ------------------------------------

{
max3_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = three_cliques, fill = graph_type), 
                   binwidth = 0.05, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = "none") +
    xlab('# of 3-cliques') 
  
max3_plot
ggsave('pictures/3cliques.png', max3_plot, limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot transitivity --------------------------------------------------------

{
trans_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = transitivity, fill = graph_type), 
                   binwidth = 0.05, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    scale_x_log10() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = "none") +
    xlab('transitivity')
  
trans_plot
ggsave('pictures/transitivity.png', trans_plot, limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot assortativity --------------------------------------------------------

{
assort_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = assortativity, fill = graph_type), 
                   binwidth = 0.02, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    theme_minimal() +
    scale_x_continuous(breaks = c(-0.4, 0, 0.5, 1)) +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = guide_legend(title = "graph type")) +
    xlab('assortativity')
  
assort_plot
ggsave('pictures/assortativity.png', assort_plot, limitsize = FALSE, 
         units = c("px"), width = 1100, height = 800, bg = 'white')
}

# plot diameter --------------------------------------------------------

{
diameter_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = diameter, fill = graph_type), 
                   binwidth = 1.2, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    scale_x_continuous() +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = "none") +
    xlab('diameter')
  
diameter_plot
ggsave('pictures/diameter.png', limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# plot meandist --------------------------------------------------------

{
meandist_plot <- all_graphs_tb %>%
    ggplot() +
    geom_histogram(aes(x = meandist, fill = graph_type), 
                   binwidth = 0.2, 
                   alpha = 0.6, color = "grey20", size = 0.05) +
    scale_x_continuous(breaks = c(5, 10, 15)) +
    theme_minimal() +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#E69F00',
                                 'grey90', '#000000'),
                      breaks = c('desmo', 'erdos', 'neuro', 
                                 'pa', 'sf')) +
    guides(fill = "none") +
    xlab('meandist')
  
meandist_plot
ggsave('pictures/meandist.png', meandist_plot, limitsize = FALSE, 
         units = c("px"), width = 800, height = 800, bg = 'white')
}

# save the data tibble ----------------------------------------------------

saveRDS(all_graphs_tb, file = "source_data/Figure4_source_data_1.RDS",
        compress = TRUE)
writeLines(capture.output(dput(all_graphs_tb)), "source_data/Figure4_source_data_1.txt")



# create multi-panel figure -----------------------------------------------
{
#combine degree and weight plot with patchwork and convert to panel with ggdraw, shared label
degree_weight <- ggdraw(degree_plot + weights_plot +
plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain')))

modul_diameter <- ggdraw(modul_plot + diameter_plot +
                          plot_layout(guides = 'collect') +
                          plot_annotation(tag_levels = list(c('C', 'D'))) & 
                          theme(plot.tag = element_text(size = 12, face='plain')))
layout = "ABCD"
mean_trans_3_ass <- ggdraw(meandist_plot + trans_plot + max3_plot + assort_plot +
                           plot_layout(design = layout, guides = 'collect') +
                           plot_annotation(tag_levels = list(c('E', 'F', 'G', 'H'))) & 
                           theme(plot.tag = element_text(size = 12, face='plain')))

layout <- "
AAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBB
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

Fig4 <- degree_weight + modul_diameter + mean_trans_3_ass +
    plot_layout(design = layout, heights = c(1, 1),
                guides = 'collect')

ggsave("figures/Figure_4.pdf", limitsize = FALSE, 
         units = c("px"), Fig4, width = 2800, height = 1300)
}

ggsave("figures/Figure_4.png", limitsize = FALSE, 
       units = c("px"), Fig4, width = 2800, height = 1300)



