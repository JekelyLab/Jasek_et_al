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