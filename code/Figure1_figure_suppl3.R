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

ggsave("figures/Figure1_figure_suppl3_desmo-tono_graph.pdf", limitsize = FALSE, 
       units = c("px"), desmo_tono_graph, width = 2000, height = 2000)

ggsave("figures/Figure1_figure_suppl3_desmo-tono_graph.png", limitsize = FALSE, 
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

tags_bf_x <- NULL
tags_bf_y <- NULL
tags_bf_z <- NULL

for (i in (1:length(treenodes_with_bf_tags_all))) {
  print(i)
  print(treenodes_with_bf_tags_all[i])
  treenodes_bf_detail <- catmaid_fetch(path = "11/treenodes/compact-detail", 
                                       body = list(treenode_ids=treenodes_with_bf_tags_all[i]))
  tags_bf_x[i] <-  treenodes_bf_detail[[1]][[3]]
  tags_bf_y[i] <-  treenodes_bf_detail[[1]][[4]]
  tags_bf_z[i] <-  treenodes_bf_detail[[1]][[5]]
}


plot_landmarks <- function() {
#plot meshes and background reference cells
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.6,
       col="grey50")
}


nopen3d(); mfrow3d(1, 2)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 600, 800)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)

#define the size of the rgl window
par3d(windowRect = c(20, 30, 1200, 800)) 

#define view and zoom
nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)


#add a text label
texts3d(35000,0, 0, text = "ventral view", col="black", cex = 2.5)

#move to next panel in rgl window
next3d(clear=F)
#define view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.55)

next3d(clear=F)

#define the lighting
clear3d(type = "lights"); 
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")
#set background color
bg3d("gray100")
next3d(clear=F)
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")

next3d(clear=F)

texts3d(35000,0, 0, text = "ventral view", col="black", cex = 2.5)

plot3d(
  desmo_x,
  desmo_y,
  desmo_z,
  add = TRUE, 
  col=hcl.colors(20000, palette='Blues'), 
  size=4, 
  alpha=0.8
)

plot_landmarks()

next3d(clear=F)

texts3d(115000,30000, 8020, text = "left view", col="black", cex = 2.5)

plot3d(
  desmo_x,
  desmo_y,
  desmo_z,
  add = TRUE, 
  col=hcl.colors(20000, palette='Blues'), 
  size=4, 
  alpha=0.8
)

plot_landmarks()

rgl.snapshot("pictures/Fig1_suppl3_desmosomes.png")

next3d(clear = T)

texts3d(35000,0, 0, text = "ventral view", col="black", cex = 2.5)

plot3d(
  tags_bf_x,
  tags_bf_y,
  tags_bf_z,
  add = TRUE, 
  col=hcl.colors(1, palette='Reds'), 
  size=5, 
  alpha=0.8
)

plot3d(
  bf_skeletons, 
  soma=T, 
  lwd=2,
  add=T,
  alpha=0.4,
  col="light pink"
) 

plot_landmarks()

next3d(clear = T)

texts3d(115000,30000, 8020, text = "left view", col="black", cex = 2.5)

plot3d(
  tags_bf_x,
  tags_bf_y,
  tags_bf_z,
  add = TRUE, 
  col=hcl.colors(1, palette='Reds'), 
  size=5, 
  alpha=0.8
)

plot3d(
  bf_skeletons, 
  soma=T, 
  lwd=2,
  add=T,
  alpha=0.4,
  col="light pink"
) 

plot_landmarks()

rgl.snapshot("pictures/Fig1_suppl3_tonofibrils.png")

next3d(clear = F)

plot3d(
  desmo_x,
  desmo_y,
  desmo_z,
  add = TRUE, 
  col=hcl.colors(20000, palette='Blues'), 
  size=4, 
  alpha=0.8
)

next3d(clear = F)

plot3d(
  desmo_x,
  desmo_y,
  desmo_z,
  add = TRUE, 
  col=hcl.colors(20000, palette='Blues'), 
  size=4, 
  alpha=0.8
)

rgl.snapshot("pictures/Fig1_suppl3_tonofibrils-desmosomes.png")

# assemble figure ---------------------------------------------------------

{
  panel_A <-desmo_tono_graph
  
  panel_B <- ggdraw() + draw_image(readPNG("pictures/Fig1_suppl3_tonofibrils.png")) 
  panel_C <- ggdraw() + draw_image(readPNG("pictures/Fig1_suppl3_desmosomes.png"))
  panel_D <- ggdraw() + draw_image(readPNG("pictures/Fig1_suppl3_tonofibrils-desmosomes.png"))
  
  layout <- "
AAAABBBB
CCCCDDDD
"
  
  Figure1_suppl3 <- panel_A + panel_B + 
    panel_C + panel_D +
    plot_layout(design = layout, heights = c(1, 1, 1, 1.6)) +
    plot_annotation(tag_levels = list(c('A', 'B',
                                        'C', 'D'))) & 
    theme(plot.tag = element_text(size = 12, face='plain'))
  
  
  ggsave("figures/Figure1_figure_supplement3.pdf", limitsize = FALSE, 
         units = c("px"), Figure1_suppl3, width = 2600, height = 2400)
  
  
  ggsave("figures/Figure1_figure_supplement3.pdf", limitsize = FALSE, 
         units = c("px"), Figure1_suppl3, width = 2600, height = 2400, bg='white')
}
