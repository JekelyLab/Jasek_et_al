skids_from_connector <- function(connector_id, pid) {
  cat_path = paste("", pid, "connectors", connector_id, sep="/")
  connectors_response = catmaid_fetch(cat_path)
  skid1 = connectors_response$partners[[1]]$skeleton_id
  connector1_type = connectors_response$partners[[1]]$relation_name
  if (length(connectors_response$partners)>1) {
    skid2 = connectors_response$partners[[2]]$skeleton_id
  } else {
    skid2 = NA
  }
  df=data.frame(connector=connector_id, x=connectors_response$x, y=connectors_response$y, z=connectors_response$z, skid1=skid1, skid2=skid2)
}



#######################################################
# desmosomes just for some groups of skids
  
muscle_neuronlist = read.neurons.catmaid("^muscle$", pid=11, 
                                           fetch.annotations = T)
  
connector_ids <- connectors(muscle_neuronlist)
  
# if some selected skids conenct with each other, those connectors will be duplicated. remove them
desmo_connector_ids <- connector_ids[connector_ids$prepost == 3,]
desmo_conenctor_ids_uniq = unique(desmo_connector_ids$connector_id)
  
t2 <- NULL
for (i in desmo_conenctor_ids_uniq) {
  print(i)
  t <- skids_from_connector(i, 11)
  t2 <- rbind(t2, t)
}

write.csv(t2, "data/desmo_mus_positions_skids.csv")

########################################################
# for all desmosomes
all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

t2 <- NULL
for (i in 1:(length(all_desmo_connectors$connectors))) {
  print(i)
  connectorid = all_desmo_connectors$connectors[[i]][[1]]
  print(connectorid)
  t <- skids_from_connector(connectorid, 11)
  t2 <- rbind(t2, t)
}

write.csv(t2, "data/desmo_all_positions_skids.csv")


########################################################
# find desmosomes connected to only one skeleton instead of 2
# there shouldn't be any of these, and any such desmosomes should be fixed in the database
desmo_with_one_partner_indices <- which(is.na(t2$connector2_type))
desmo_with_one_partner_ids <- desmo_conenctor_ids_uniq[desmo_with_one_partner_indices]
