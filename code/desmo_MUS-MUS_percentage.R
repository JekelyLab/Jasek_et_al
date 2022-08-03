library(natverse)

source("~/R/conn.R")

all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

desmo_count_MUS = 0
desmo_count_other = 0

for (i in (1:length(all_desmo_connectors$connectors))) {
  print(i)
  skid1 <- all_desmo_connectors$partners[[i]][[1]][[3]]
  skid2 <- all_desmo_connectors$partners[[i]][[2]][[3]]
  skid1annotations <- catmaid_get_annotations_for_skeletons(skid1, pid=11)
  skid2annotations <- catmaid_get_annotations_for_skeletons(skid2, pid=11)
  if ("muscle" %in% skid1annotations$annotation & "muscle" %in% skid2annotations$annotation ) {
    desmo_count_MUS = desmo_count_MUS + 1
    print(paste("skids ", skid1, " and ", skid2))
  }
  else {
    desmo_count_other = desmo_count_other + 1
  }
}

percentage <- desmo_count_MUS / (desmo_count_other + desmo_count_MUS)
