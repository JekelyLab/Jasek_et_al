#  don't source the giant packages and connection script, it's too big and almost all of it is not needed

library(catmaid)
library(tidyverse)

source("~/R/conn.R")

################################################################################
# numbers of cells with tonofibrils, desmosomes or both, by cell type


# get all skids with "black fibers" tag
# there are also annotations for this but tags are more reliable because they
# will always stay with the structure when splitting neurons

tags_all <- catmaid_get_label_stats(pid = 11)

skids_with_bf_tags_all <- tags_all %>% 
  filter(labelName=="black fibers") %>% 
  select(skeletonID) %>%
  unlist() %>%
  unique()


# get all skids with desmosomes
desmo_connectors_all <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

partners1 <- as.data.frame(sapply(desmo_connectors_all$partners, "[[", 1))
partners1 <- as.vector(sapply(partners1, "[[", 3))

partners2 <- as.data.frame(sapply(desmo_connectors_all$partners, "[[", 2))
partners2 <- as.vector(sapply(partners2, "[[", 3))

skids_desmo_all <- unique(c(partners1, partners2))


# function to count number of cells with desmosomes, tonofibrils, both or neither, per celltype
count_desmo_tono <- function(annot) {
  annotation_query = paste("^", annot,"$", sep="")
  cells_with_annotation <<- catmaid_query_by_annotation(annotation_query, 
                                                       pid = 11)
  
  skids_with_annotation <<- cells_with_annotation$skid
  
  skids_desmo <<- intersect(skids_with_annotation, skids_desmo_all)
  skids_nondesmo <<- setdiff(skids_with_annotation, skids_desmo)
  skids_bf <<- intersect(skids_with_annotation, skids_with_bf_tags_all)
  skids_desmo_bf <<- intersect(skids_desmo, skids_bf)
  
  skids_nonbf <<- setdiff(skids_with_annotation, skids_bf)
  skids_desmo_nonbf <<- intersect(skids_desmo, skids_nonbf)
  
  skids_nondesmo_bf <<- intersect(skids_nondesmo, skids_bf)
  
  skids_nondesmo_nonbf <<- intersect(skids_nondesmo, skids_nonbf)
 
}



celltypes_bf_desmo <- data.frame()

# non-neuronal cell types excluding somatic muscles
for (i in c(1:36, 79, 90:92)){
  annotation = paste("celltype_non_neuronal", i, sep="")
  count_desmo_tono(annotation)
  
  df <- data.frame(celltype=annotation,
                   desmo=length(skids_desmo_nonbf),
                   desmo_tonofibrils=length(skids_desmo_bf),
                   tonofibrils=length(skids_nondesmo_bf),
                   neither=length(skids_nondesmo_nonbf),
                   total=length(skids_with_annotation))
  celltypes_bf_desmo <- rbind(celltypes_bf_desmo, df)
}


celltype_names <- read.csv("data/non_neuronal_celltypes_names.csv")

celltypes_bf_desmo_with_names <- right_join(celltype_names,
                                            celltypes_bf_desmo,
                                            by = c("CATMAID.annotation" = "celltype"))


# somatic muscles
len_skids_desmo_nonbf <- len_skids_desmo_bf <- len_skids_nondesmo_bf <- len_skids_nondesmo_nonbf <- len_skids_with_annotation <- 0
for (i in c(37:78, 80:89)) {
  annotation = paste("celltype_non_neuronal", i, sep="")
  count_desmo_tono(annotation)
  
  len_skids_desmo_nonbf <- len_skids_desmo_nonbf + length(skids_desmo_nonbf)
  len_skids_desmo_bf <- len_skids_desmo_bf + length(skids_desmo_bf)
  len_skids_nondesmo_bf <- len_skids_nondesmo_bf + length(skids_nondesmo_bf)
  len_skids_nondesmo_nonbf <- len_skids_nondesmo_nonbf + length(skids_nondesmo_nonbf)
  len_skids_with_annotation <- len_skids_with_annotation + length(skids_with_annotation)
  
}
df=data.frame(Name="somatic muscles",
              CATMAID.annotation="celltype_non_neuronal37-78,80-89",
              desmo=len_skids_desmo_nonbf,
              desmo_tonofibrils=len_skids_desmo_bf,
              tonofibrils=len_skids_nondesmo_bf,
              neither=len_skids_nondesmo_nonbf,
              total=len_skids_with_annotation)
celltypes_bf_desmo_with_names <- rbind(celltypes_bf_desmo_with_names, df)

# neurons
len_skids_desmo_nonbf <- len_skids_desmo_bf <- len_skids_nondesmo_bf <- len_skids_nondesmo_nonbf <- len_skids_with_annotation <- 0
for (i in c(1:200)) {
  annotation = paste("celltype", i, sep="")
  count_desmo_tono(annotation)
  
  len_skids_desmo_nonbf <- len_skids_desmo_nonbf + length(skids_desmo_nonbf)
  len_skids_desmo_bf <- len_skids_desmo_bf + length(skids_desmo_bf)
  len_skids_nondesmo_bf <- len_skids_nondesmo_bf + length(skids_nondesmo_bf)
  len_skids_nondesmo_nonbf <- len_skids_nondesmo_nonbf + length(skids_nondesmo_nonbf)
  len_skids_with_annotation <- len_skids_with_annotation + length(skids_with_annotation)
  
}
df=data.frame(Name="neurons",
              CATMAID.annotation="celltype1-200",
              desmo=len_skids_desmo_nonbf,
              desmo_tonofibrils=len_skids_desmo_bf,
              tonofibrils=len_skids_nondesmo_bf,
              neither=len_skids_nondesmo_nonbf,
              total=len_skids_with_annotation)
celltypes_bf_desmo_with_names <- rbind(celltypes_bf_desmo_with_names, df)


write.csv(celltypes_bf_desmo_with_names, "data/percent_cells_with_desmo_bf_by_celltype.csv", row.names = FALSE)

celltypes_bf_desmo_with_names_arranged <- arrange(celltypes_bf_desmo_with_names,
                                                  neither/total,
                                                  desc(desmo/total),
                                                  desc(desmo_tonofibrils/total),
                                                  desc(tonofibrils/total))

celltypes_bf_desmo_with_names_tidy <-  celltypes_bf_desmo_with_names_arranged %>%
  select(-total) %>%
  pivot_longer(
    cols = c("desmo", "desmo_tonofibrils", "tonofibrils", "neither"), 
    names_to = "characteristic", 
    values_to = "count")

desmo_tono_graph <- celltypes_bf_desmo_with_names_tidy %>%
  ggplot(aes(Name,
             count,
             fill=factor(characteristic,
                         levels=c("desmo", "desmo_tonofibrils", "tonofibrils", "neither")))) +
  geom_bar(position="fill",
           stat = "identity") +
  scale_x_discrete(limits = rev(celltypes_bf_desmo_with_names_arranged$Name)) +
  scale_y_reverse() +
  geom_text(aes(label = count),
            position = "fill",
            hjust = 1,
            size = 2,
            family="Arial") +
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
        text=element_text(family="Arial"))

ggsave("figures/Figure1_figure_suppl3.pdf", limitsize = FALSE, 
       units = c("px"), desmo_tono_graph, width = 2000, height = 2000)

ggsave("figures/Figure1_figure_suppl3.png", limitsize = FALSE, 
       units = c("px"), desmo_tono_graph, width = 2000, height = 2000, bg = 'white')



################################################################################
# visualize desmosomes, "black fibers" tags, and skeletons with "black fibers" annotations in 3D

desmo_connectors_all <- catmaid_fetch(path = "11/connectors/", body = list(relation_type="desmosome_with", with_partners="true"))

desmo_x <- sapply(desmo_connectors_all$connectors, "[[", 2)
desmo_y <- sapply(desmo_connectors_all$connectors, "[[", 3)
desmo_z <- sapply(desmo_connectors_all$connectors, "[[", 4)

tags_all <- catmaid_get_label_stats(pid = 11)

skids_with_bf_tags_all <- tags_all %>% 
  filter(labelName=="black fibers") %>% 
  select(skeletonID) %>%
  unlist() %>%
  unique()

treenodes_with_bf_tags_all <- tags_all %>% 
  filter(labelName=="black fibers") %>% 
  select(treenodeID) %>%
  unlist() %>%
  unique()


bf_skeletons = nlapply(
  read.neurons.catmaid(
    "^black fibers$", 
    pid=11, 
    fetch.annotations = T
  ), 
  function(x) smooth_neuron(x, sigma=6000)
)

nopen3d() # opens a pannable 3d window
mfrow3d(1, 1)  #defines the two scenes
par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
par3d(zoom=0.52)

plot3d(
  outline, 
  add=T, 
  alpha=0.04,
  col="#E2E2E2"
) 

plot3d(
  yolk, 
  add=T, 
  alpha=0.06,
  col="#E2E2E2"
) 

plot3d(
  acicula, 
  soma=T, 
  lwd=3,
  add=T, 
  alpha=0.5,
  col="black"
) 

plot3d(
  desmo_x,
  desmo_y,
  desmo_z,
  add = TRUE, 
  col=hcl.colors(20000, palette='Oranges'), 
  size=3, 
  alpha=0.5
)

plot3d(
  bf_skeletons, 
  soma=T, 
  lwd=2,
  add=T,
  alpha=0.3,
  col="light blue"
) 

for (i in c(treenode_ids=treenodes_with_bf_tags_all)) {
  treenodes_bf_detail <- catmaid_fetch(path = "11/treenodes/compact-detail", body = list(treenode_ids=i))
  
  plot3d(
    treenodes_bf_detail[[1]][[3]],
    treenodes_bf_detail[[1]][[4]],
    treenodes_bf_detail[[1]][[5]],
    add = TRUE, 
    col=hcl.colors(1, palette='Blues'), 
    size=3, 
    alpha=0.5
  )
}

