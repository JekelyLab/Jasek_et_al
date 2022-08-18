library(natverse)

source("~/R/conn.R")

all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

desmo_count_MUS = 0
desmo_count_MUS_same = 0
desmo_count_MUS_other = 0
desmo_count_other = 0
desmo_count_lamina = 0

for (i in (1:length(all_desmo_connectors$connectors))) {
  print(i)
  skid1 <- all_desmo_connectors$partners[[i]][[1]][[3]]
  skid2 <- all_desmo_connectors$partners[[i]][[2]][[3]]
  skid1annotations <- catmaid_get_annotations_for_skeletons(skid1, pid=11)
  skid2annotations <- catmaid_get_annotations_for_skeletons(skid2, pid=11)
  mus1 <- skid1annotations$annotation[grepl("^muscle$", skid1annotations$annotation)]
  mus2 <- skid2annotations$annotation[grepl("^muscle$", skid2annotations$annotation)]
  lamina1 <- skid1annotations$annotation[grepl("^basal lamina$", skid1annotations$annotation)]
  lamina2 <- skid2annotations$annotation[grepl("^basal lamina$", skid2annotations$annotation)]
  if (length(mus1)>0 & length(mus2)>0) {
    desmo_count_MUS = desmo_count_MUS + 1
    print(paste("skids ", skid1, " and ", skid2, "are muscles"))
    celltype <- skid1annotations$annotation[grepl("celltype_non_neuronal.", skid1annotations$annotation)]
    print(celltype)
    if (length(celltype)>0) {
      if (celltype %in% skid2annotations$annotation) {
        desmo_count_MUS_same = desmo_count_MUS_same + 1
        print("same muscle type")
      }
      else {
        desmo_count_MUS_other = desmo_count_MUS_other + 1
        print("different muscle types")
      }
    }
    else {
      desmo_count_MUS_other = desmo_count_MUS_other + 1
    }
  }
  else if (length(lamina1)>0 | length(lamina2)>0) { 
    print(paste("skids", skid1, "or", skid2, "are basal lamina"))
    desmo_count_lamina = desmo_count_lamina + 1
  }
  else {
    desmo_count_other = desmo_count_other + 1
  }
}


percentage_hemidesmo <- desmo_count_lamina / length(all_desmo_connectors$connectors)

# percentage of desmosoems which are muscle to muscle
percentage_MUS_desmo <- desmo_count_MUS / length(all_desmo_connectors$connectors)

# percentage of muscle to muscle desmosoems which are between same muscle types
percentage_MUS_same_desmo <- desmo_count_MUS_same / (desmo_count_MUS_other + desmo_count_MUS_same)

# percentage of muscle to muscle desmosoems which are between different muscle types
percentage_MUS_other_desmo <- desmo_count_MUS_other / (desmo_count_MUS_other + desmo_count_MUS_same)
