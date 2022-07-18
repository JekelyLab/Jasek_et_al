skids_from_connector <- function(connector_id, pid) {
  # TODO: make this work for synapses with multiple post synaptic partners
  # but I don't care about synapses right now
  cat_path = paste("", pid, "connectors", connector_id, sep="/")
  connectors_response = catmaid_fetch(cat_path)
  skid1 = connectors_response$partners[[1]]$skeleton_id
  connector1_type = connectors_response$partners[[1]]$relation_name
  if (length(connectors_response$partners)>1) {
    skid2 = connectors_response$partners[[2]]$skeleton_id
    connector2_type = connectors_response$partners[[2]]$relation_name
  } else {
    skid2 = NA
    connector2_type = NA
  }
  df=data.frame(skid1, connector1_type, skid2, connector2_type)
}


MUSlongD = nlapply(read.neurons.catmaid("^celltype_non_neuronal75$", pid=11, 
                                      fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
#get the connectors
MUSconn <- connectors(MUSlongD)

MUSconn3[207,2]
skids_from_connector(7047057, 11)
skids_from_connector(7003151, 11)

#separate synapse (prepost 1) and desmosomal (prepost 3) connectors
MUSconn1 <- MUSconn[MUSconn$prepost %in% c("1"),]
MUSconn3 <- MUSconn[MUSconn$prepost %in% c("3"),]
(MUSconn3[1:1000,2])
MUSconn3[208,2]

#iterate through the connector list and run function
sapply(MUSconn3[,2], function(x) skids_from_connector(x, 11))

(MUSconn3[546,2])

desmo_list <- list()
for (i in 1:10) {
  print(i)
  desmo_list[[i]] <- skids_from_connector(MUSconn3[i,2], 11)
}
desmo_list[[1]][1]
desmo_list[[1]][3]
desmo_list
library(igraph)
g <- make_empty_graph() %>%
add_vertices(1, attr = list(name = 1853803)) %>%
add_vertices(1, attr = list(name = 137882)) 
V(g)
E(g)

library(tidygraph)
g %>%
  as_tbl_graph()%>% 
  bind_edges(data.frame(from = '1853803', to = '137882'), node_key = 'name')


