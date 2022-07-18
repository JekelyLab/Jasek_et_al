#R code to generate the graph in Figure2 panel I on the segmentaldistributionof muscles in Jasek et al. 2021 Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely

source("code/Packages_and_Connection.R")


#############################
#define function to retrieve skids from a neuron list based on one to three annotations
skids_by_annotation <- function(neuron_list,annotation1,annotation2,annotation3){
  skids1 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation1,1]))
  if(missing(annotation2)){return(skids1) #if annotation2 is missing, will return skids matching the first annotation
  } else {
    skids2 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation2,1]))
  }
  if(missing(annotation3)){return(intersect(skids1,skids2)) #if annotation3 is missing, will return skids matching annotations 1 and 2
  }   else  {
    skids3 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation3,1]))
  }
  skids1_2 <- intersect(skids1,skids2)
  return(intersect(skids1_2,skids3)) #will return the shared skids between the three annotations  
}



annotation_non_neuronal_celltypelist = list()
#we read all non-neuronal celltypes - muscles from 37-89 and all annotations
for (i in c(37:89)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_non_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define the six body regions, matching the catmaid annotations
regions <- c('episphere','segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')
sides <- c('left_side','right_side')
types <- paste('celltype_non_neuronal', 37:89, sep='')

#we retrieve those skids that match our annotations
#three iterated lapply functions to retrieve skids by muscle type, body region and body side
muscle_per_side_per_body_region <- lapply(sides, function(s) lapply(regions, function(r) lapply(types, function(m) skids_by_annotation(annotation_non_neuronal_celltypelist,s,m,r))))

n_cells_left <- matrix(nrow=6,ncol=53);n_cells_right <- matrix(nrow=6,ncol=53)
#count the occurrence of each type in each body region per side
for (j in 1:6){for (i in 1:53){
n_cells_left[j,i] <- length(muscle_per_side_per_body_region[[1]][[j]][[i]])
n_cells_right[j,i] <- length(muscle_per_side_per_body_region[[2]][[j]][[i]])
}}

library(heatmaply)
#add row and column names
row.names(n_cells_left) <- c('head left','sg0 left','sg1 left','sg2 left','sg3 left','pyg left')
row.names(n_cells_right) <- c('head right','sg0 right','sg1 right','sg2 right','sg3 right','pyg right')
#combine left and right matrix
n_cells <- rbind(n_cells_left,n_cells_right)

#celltype name lists
non_neuronal_celltype_names = c("akrotroch", "crescent cell", "prototroch", "nuchal cilia", "metatroch", "paratroch", "spinGland", "covercell", "ciliatedGland", "eyespot pigment cell", "pigment cell AE", "bright droplets parapodial", "macrophage", "yolk cover cell", "flat glia", "EC rad glia ", "MVGland", "microvillarCell", "protonephridium", "nephridium", "nephridiumTip", "chaeta", "aciculoblast", "circumacicular", "hemichaetal", "ER circumchaetal", "noER circumchaetal", "EC circumchaetal", "HeadGland", "InterparaGland", "spinMicroGland", "CB pigment", "vacuolar cell_head", "Glia pigmented", "pygidial pigment", "meso", "MUSac_notA", "MUSac_notP", "MUSac_notM", "MUSac_neuAV", "MUSac_neuPD", "MUSac_neuPV", "MUSac_neuDy", "MUSac_neuDx", "MUSac_neuDach", "MUSac_neure", "MUSac_i", "MUSob-ant_re", "MUSob-ant_arc", "MUSob-ant_m-pp", "MUSob-ant_ml-pp", "MUSob-ant_l-pp", "MUSob-ant_trans", "MUSob-post_notD", "MUSob-post_neuDlong", "MUSob-post_neuDprox", "MUSob-post_neuDdist", "MUSob-post_neuV", "MUSob-post_notV", "MUSob-postM", "MUSob-post_noty", "MUSob-post_i", "MUSchae_notDob", "MUSchae_notD", "MUSchae_notDn", "MUSchae_notA", "MUSchae_notAac", "MUSchae_notAre", "MUSchae_neuVob", "MUSchae_neuDac", "MUSchae_neuAVo", "MUSchae_neuAVt", "MUSchae_Are", "MUStrans_pyg", "MUSlong_D", "MUSlong_V", "MUSax", "MUSring", "MUSph", "MUSll", "MUSant", "MUSly", "MUSpl", "MUSci", "MUSch", "MUSpx", "MUSpr", "MUStri", "MUSmed_head", "EC")
#these are the muscle cell type names
colnames(n_cells)<-non_neuronal_celltype_names[37:89]

#initialise new active plotting device
plot.new()
#plot heatmap
hm <- heatmaply(n_cells,Rowv=F, Colv=F,
          cellnote=ifelse(n_cells==0, NA, n_cells),
          cellnote_size=17,
          cellnote_textposition='middle center',
          col=c('white','grey20','#1098CD','red'),
          column_text_angle = 90,
          fontsize_col = 12,fontsize_row = 10)
hm
#save table
write.csv(n_cells, file = "source_data/Figure2I_source_data_1.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


#save output of viewer to html using the Export/save as web page drop down menu in the Viewer

#save as html
webshot2::webshot(url="pictures/Muscle_Celltypes_numbers.html",
                  file="pictures/Muscle_Celltypes_numbers.png",
                  cliprect = c(94,0,1415, 530), zoom=5,
                  vwidth = 1600,
                  vheight = 744)

