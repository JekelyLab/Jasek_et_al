#  don't source the giant packages and connection script, it's too big and almost all of it is not needed

library(catmaid)
library(tidyverse)

source("~/R/conn.R")

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

celltypes_bf_desmo_with_names_arranged <- arrange(celltypes_bf_desmo_with_names,
                                                  other/total,
                                                  desc(desmo/total),
                                                  desc(desmo_tonofibrils/total),
                                                  desc(tonofibrils/total))

celltypes_bf_desmo_with_names_tidy <-  celltypes_bf_desmo_with_names_arranged %>%
  select(-total) %>%
  pivot_longer(
    cols = c("desmo", "desmo_tonofibrils", "tonofibrils", "other"), 
    names_to = "characteristic", 
    values_to = "count")

desmo_tono_graph <- celltypes_bf_desmo_with_names_tidy %>%
  ggplot(aes(Name, count, fill=factor(characteristic, levels=c("desmo", "desmo_tonofibrils", "tonofibrils", "other")))) +
  geom_bar(position="fill", stat = "identity") +
  scale_x_discrete(limits = rev(celltypes_bf_desmo_with_names_arranged$Name)) +
  scale_y_reverse() +
  geom_text(aes(label = count), position = "fill", hjust=1, size = 2, family="ArialMT") +
  coord_flip() + 
  #scale_fill_manual(values = c("#E69F00","#D55E00","#CC79A7","lightgrey")) + # Okabe Ito
  scale_fill_manual(values = c("#4477AA","#AA3377","#EE6677","lightgrey")) + # Tol
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        text=element_text(family="ArialMT"))

ggsave("figures/Figure1_figure_suppl3.pdf", limitsize = FALSE, 
       units = c("px"), desmo_tono_graph, width = 2000, height = 2000)

ggsave("figures/Figure1_figure_suppl3.png", limitsize = FALSE, 
       units = c("px"), desmo_tono_graph, width = 2000, height = 2000, bg = 'white')
