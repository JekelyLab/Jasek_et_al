library(catmaid)
library(tidyverse)

source("~/R/conn.R")

################################################################################
# Calculate percentages of desmosomes which are:
#     - hemidesmosomes
#     - muscle to muscle
#     - muscle to muscle, same muscle type
#     - muscle to muscle, different muscle types


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



################################################################################
# Percentage of cells with desmosomes, which also have black fibers (intermediate filaments)

# get info about all desmo connectors in the Naomi project
all_desmo_connectors <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

# extract ids of desmo connectors
all_desmo_ids <- sapply(all_desmo_connectors$connectors, "[[", 1)

# find skids connected to the desmosomes
partners1 <- as.data.frame(sapply(all_desmo_connectors$partners, "[[", 1))
partners1 <- as.vector(sapply(partners1, "[[", 3))

partners2 <- as.data.frame(sapply(all_desmo_connectors$partners, "[[", 1))
partners2 <- as.vector(sapply(partners2, "[[", 3))

all_skids <- unique(c(partners1, partners2))

annotations_all = catmaid_get_annotations_for_skeletons(all_skids, pid=11)

annotation_mus_skids = annotations_all[grepl("^muscle$", annotations_all$annotation),1]
annotation_nonmus_skids = setdiff(all_skids, annotation_mus_skids)
annotation_bl_skids = annotations_all[grepl("^basal lamina$", annotations_all$annotation),1]
annotation_bf_skids = annotations_all[grepl("^black fibers$", annotations_all$annotation),1]
annotation_bf_nonbl_skids = setdiff(annotation_bf_skids, annotation_bl_skids)

num_skids_mus_bf = length(intersect(annotation_mus_skids, annotation_bf_skids))
num_skids_nonmus_bf = length(intersect(annotation_nonmus_skids, annotation_bf_skids))
num_nonmus = length(annotation_nonmus_skids)

annotation_nonmus_nonbf = setdiff(annotation_nonmus_skids, annotation_bf_skids)


################################################################################
# numbers of cells with tonofibrils, desmosomes or both, by cell type

celltypes_bf_desmo <- data.frame()

# non-neuronal cell types excluding somatic muscles
for (i in c(1:36, 79, 90:92)){
#for (i in c(1:92)){
  #annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  annotation = paste("celltype_non_neuronal", i, sep="")
  neurons <- read.neurons.catmaid(annotation, pid=11) # I'm not fetching annotations here, but doing it in a separate step later to limit the size of requests. Annotations with over a thousand skids can cause problems
  skids <- as.integer(sapply(neurons, "[[", 1))
  connectors <- connectors(neurons)
  if (length(connectors$skid) == 0) {
    skids_desmo <- ""
  } else {
    desmo <- connectors[connectors$prepost %in% 3,]
    skids_desmo <- unique(desmo$skid)
  }
  annotations <- catmaid_get_annotations_for_skeletons(skids, pid=11)
  #annotations <- attr(neurons, 'anndf')
  skids_desmo <- as.integer(skids_desmo)
  skids_nondesmo <- setdiff(skids, skids_desmo)
  skids_bf <- annotations[grepl("^black fibers$", annotations$annotation), "skid"]
  skids_desmo_bf <- intersect(skids_desmo, skids_bf)
  
  skids_nonbf <- setdiff(skids, skids_bf)
  skids_desmo_nonbf <- intersect(skids_desmo, skids_nonbf)
  
  skids_nondesmo_bf <- intersect(skids_nondesmo, skids_bf)
  
  skids_nondesmo_nonbf <- intersect(skids_nondesmo, skids_nonbf)
  
  df=data.frame(celltype=annotation, desmo=length(skids_desmo_nonbf), desmo_tonofibrils=length(skids_desmo_bf), tonofibrils=length(skids_nondesmo_bf), other=length(skids_nondesmo_nonbf), total=length(skids))
  celltypes_bf_desmo <- rbind(celltypes_bf_desmo, df)
}

celltype_names <- read.csv("data/non_neuronal_celltypes_names.csv")
 
celltypes_bf_desmo_with_names <- right_join(celltype_names, celltypes_bf_desmo, by = c("CATMAID.annotation" = "celltype"))
#celltypes_bf_desmo <- celltypes_bf_desmo %>% remove_rownames %>% column_to_rownames(var="CATMAID.annotation")


# somatic muscles
len_skids_desmo_nonbf <- len_skids_desmo_bf <- len_skids_nondesmo_bf <- len_skids_nondesmo_nonbf <- len_skids <- 0
for (i in c(37:78, 80:89)) {
  annotation = paste("celltype_non_neuronal", i, sep="")
  neurons <- read.neurons.catmaid(annotation, pid=11) # I'm not fetching annotations here, but doing it in a separate step later to limit the size of requests. Annotations with over a thousand skids can cause problems
  skids <- as.integer(sapply(neurons, "[[", 1))
  connectors <- connectors(neurons)
  if (length(connectors$skid) == 0) {
    skids_desmo <- ""
  } else {
    desmo <- connectors[connectors$prepost %in% 3,]
    skids_desmo <- unique(desmo$skid)
  }
  annotations <- catmaid_get_annotations_for_skeletons(skids, pid=11)
  skids_desmo <- as.integer(skids_desmo)
  skids_nondesmo <- setdiff(skids, skids_desmo)
  skids_bf <- annotations[grepl("^black fibers$", annotations$annotation), "skid"]
  skids_desmo_bf <- intersect(skids_desmo, skids_bf)
  
  skids_nonbf <- setdiff(skids, skids_bf)
  skids_desmo_nonbf <- intersect(skids_desmo, skids_nonbf)
  
  skids_nondesmo_bf <- intersect(skids_nondesmo, skids_bf)
  
  skids_nondesmo_nonbf <- intersect(skids_nondesmo, skids_nonbf)
  
  len_skids_desmo_nonbf <- len_skids_desmo_nonbf + length(skids_desmo_nonbf)
  len_skids_desmo_bf <- len_skids_desmo_bf + length(skids_desmo_bf)
  len_skids_nondesmo_bf <- len_skids_nondesmo_bf + length(skids_nondesmo_bf)
  len_skids_nondesmo_nonbf <- len_skids_nondesmo_nonbf + length(skids_nondesmo_nonbf)
  len_skids <- len_skids + length(skids)
  
}
df=data.frame(Name="somatic muscles", CATMAID.annotation="celltype_non_neuronal37-78,80-89", desmo=len_skids_desmo_nonbf, desmo_tonofibrils=len_skids_desmo_bf, tonofibrils=len_skids_nondesmo_bf, other=len_skids_nondesmo_nonbf, total=len_skids)
celltypes_bf_desmo_with_names <- rbind(celltypes_bf_desmo_with_names, df)

# neurons
len_skids_desmo_nonbf <- len_skids_desmo_bf <- len_skids_nondesmo_bf <- len_skids_nondesmo_nonbf <- len_skids <- 0
for (i in c(1:200)) {
  annotation = paste("celltype", i, sep="")
  neurons <- read.neurons.catmaid(annotation, pid=11) # I'm not fetching annotations here, but doing it in a separate step later to limit the size of requests. Annotations with over a thousand skids can cause problems
  skids <- as.integer(sapply(neurons, "[[", 1))
  connectors <- connectors(neurons)
  if (length(connectors$skid) == 0) {
    skids_desmo <- ""
  } else {
    desmo <- connectors[connectors$prepost %in% 3,]
    skids_desmo <- unique(desmo$skid)
  }
  annotations <- catmaid_get_annotations_for_skeletons(skids, pid=11)
  skids_desmo <- as.integer(skids_desmo)
  skids_nondesmo <- setdiff(skids, skids_desmo)
  skids_bf <- annotations[grepl("^black fibers$", annotations$annotation), "skid"]
  skids_desmo_bf <- intersect(skids_desmo, skids_bf)
  
  skids_nonbf <- setdiff(skids, skids_bf)
  skids_desmo_nonbf <- intersect(skids_desmo, skids_nonbf)
  
  skids_nondesmo_bf <- intersect(skids_nondesmo, skids_bf)
  
  skids_nondesmo_nonbf <- intersect(skids_nondesmo, skids_nonbf)
  
  len_skids_desmo_nonbf <- len_skids_desmo_nonbf + length(skids_desmo_nonbf)
  len_skids_desmo_bf <- len_skids_desmo_bf + length(skids_desmo_bf)
  len_skids_nondesmo_bf <- len_skids_nondesmo_bf + length(skids_nondesmo_bf)
  len_skids_nondesmo_nonbf <- len_skids_nondesmo_nonbf + length(skids_nondesmo_nonbf)
  len_skids <- len_skids + length(skids)
  
}
df=data.frame(Name="neurons", CATMAID.annotation="celltype1-200", desmo=len_skids_desmo_nonbf, desmo_tonofibrils=len_skids_desmo_bf, tonofibrils=len_skids_nondesmo_bf, other=len_skids_nondesmo_nonbf, total=len_skids)
celltypes_bf_desmo_with_names <- rbind(celltypes_bf_desmo_with_names, df)


write.csv(celltypes_bf_desmo_with_names, "data/percent_cells_with_desmo_bf_by_celltype.csv", row.names = FALSE)

