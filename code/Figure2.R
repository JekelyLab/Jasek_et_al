#Desmosomal connectome paper R - Natverse code to plot all muscle cell types
#code for Figure2 of the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely

source("code/Packages_and_Connection.R")

# read skeletons ----------------------------------------------------------

{
#read some meshes and cells from catmaid for plotting an anatomical background
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                              invertFaces = T, conn = NULL, pid = 11)

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)

acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11, 
                                             fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

chaeta = nlapply(read.neurons.catmaid("^celltype_non_neuronal22$", pid=11, 
                                              fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

}

#read muscles by cell type and body side

{
celltypelist = list()
annotation_celltypelist = list()

#read all non-neuronal celltypes from 37-89 (the muscles) and all annotations
for (i in c(37:89)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation
  celltypelist[[i]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = T) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

left_skids = list(); counter=0
for (df1 in annotation_celltypelist){    #iterate through the cell type list 
  counter=counter+1
  #add cells with left_side annotation to left_skids list
  left_skids[[counter]] <- df1[df1$annotation == "left_side",1]
}

}

# 3d plotting  ------------------------------------------------------------

#define randomised color list for muscle groups
{
cols_obant = hcl.colors(16, palette='Purple-Blue')
color_obant = sample(cols_obant[1:8])
cols_obpost = hcl.colors(16, palette='Green-Yellow')
color_obpost = sample(cols_obpost[1:9])
cols_ac = hcl.colors(20, palette='Blues 3')
color_ac = sample(cols_ac[1:14])
cols_ch = hcl.colors(16, palette='OrRd')
color_ch = sample(cols_ch[1:14])
cols_tr = hcl.colors(10, palette='YlOrBr')
color_tr = sample(cols_tr[1:1])
cols_lon = hcl.colors(10, palette='Terrain 2')
color_lon = sample(cols_lon[1:3])
cols_head = hcl.colors(20, palette='Warm')
color_head = sample(cols_head[1:12])
}

plot_two_panel_background <- function() {
# 3d plotting of cells
nopen3d(); mfrow3d(1, 2)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 600, 800)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)

#define the size of the rgl window
par3d(windowRect = c(20, 30, 1200, 800)) 

#define view and zoom
nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)

#plot meshes and background reference cells
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey50")

#chatetae were only plotted for the chaetal muscle group
#plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.4,col="black")

#add a text label
texts3d(35000,0, 0, text = "ventral view", col="black", cex = 2.5)

#move to next panel in rgl window
next3d(clear=F)
#define view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.55)

#plot meshes and background reference cells to the second panel
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey50")

#optional
#plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.4,col="black")

#add a text label
texts3d(115000,30000, 8020, text = "left view", col="black", cex = 2.5)
next3d(clear=F)

#define the lighting
clear3d(type = "lights"); 
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")
#set background color
bg3d("gray100")
next3d(clear=F)
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")
}

#these are the annotation numbers for the different muscle celltypes in catmaid, same in our lists
#acicular 37:47 - Blues 3
#chaetal 63:73 - OrRd
#oblique anterior 48:53 - Purple-Blue
#oblique posterior 54:62 - Emrld
#transverse 74 - brown
#longitudinal 75-77 - Terrain 2
#digestive-head 78-89 - Warm

#plot all cells
plot_two_panel_background()

{
counter = 0
for (j in c(37:89)){    #iterate through the cell type list 
  next3d(clear=F); counter = counter+1; print(j)
  nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)
  
  #assign color
  if (j>36 & j<48){color=color_ac}
  else if (j>62 & j<74){color=color_ch}
  else if (j>47 & j<54){color=color_obant}
  else if (j>53 & j<63){color=color_obpost}
  else if (j==74){color= color_tr}
  else if (j>74 & j<78){color= color_lon}
  else if (j>77 & j<90){color= color_head}
  
  if (j==37 || j==48 || j==54 || j==63 || j==74 || j==75 || j==78 ){counter=1}
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")

  next3d(clear=F)
  nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1)); par3d(zoom=0.55)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle cell types are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}
#make snapshot
rgl.snapshot("pictures/MUScelltypes_all.png")

close3d()
}

#plot acicular muscles
{
plot_two_panel_background()
next3d(clear=F)
plot3d(scalebar_50um_ventral, color = 'black', lwd = 2)
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
counter = 0
for (j in c(37:47)){    #iterate through the celltype list 
  counter = counter+1; print(j)
   #assign color
  color= color_ac
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  next3d(clear=F)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}

rgl.snapshot("pictures/MUScelltypes_ac.png")
close3d()
}

#plot anterior oblique muscles
{
plot_two_panel_background()
next3d(clear=F)
counter = 0
for (j in c(48:53)){    #iterate through the celltype list 
    counter = counter+1; print(j)
    #assign color
    color= color_obant
    cells=celltypelist[[j]]
    #smooth with sigma 6000
    cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
    #plot cells
    plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
    next3d(clear=F)
    skids_raw=left_skids[j]
    #retrieve from catmaid all left cells of a type
    for (skids in skids_raw){
      if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
      cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                           fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
      plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
    }
  }
  
  rgl.snapshot("pictures/MUScelltypes_ant_ob.png")
  close3d()
}

#plot posterior oblique muscles
{
plot_two_panel_background()
next3d(clear=F)
counter = 0
for (j in c(54:62)){    #iterate through the celltype list 
  counter = counter+1; print(j)
  #assign color
  color= color_obpost
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  next3d(clear=F)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}

rgl.snapshot("pictures/MUScelltypes_post_ob.png")
close3d()
}

#plot chaetal muscles
{
plot_two_panel_background()
  next3d(clear=F)
  counter = 0
  for (j in c(63:73)){    #iterate through the celltype list 
    counter = counter+1; print(j)
    #assign color
    color= color_ch
    cells=celltypelist[[j]]
    #smooth with sigma 6000
    cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
    #plot cells
    plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
    next3d(clear=F)
    skids_raw=left_skids[j]
    #retrieve from catmaid all left cells of a type
    for (skids in skids_raw){
      if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
      cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                           fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
      plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
    }
  }
  
  rgl.snapshot("pictures/MUScelltypes_chae.png")
  close3d()
}


#plot transverse muscles
{
plot_two_panel_background()
next3d(clear=F)
counter = 0
for (j in c(74)){    #iterate through the celltype list 
  counter = counter+1; print(j)
  #assign color
  color= color_tr
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  next3d(clear=F)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}

rgl.snapshot("pictures/MUScelltypes_trans.png")
close3d()
}

#plot Longitudinal muscles
{
plot_two_panel_background()
next3d(clear=F)
counter = 0
for (j in c(75:77)){    #iterate through the celltype list 
  counter = counter+1; print(j)
  #assign color
  color= color_lon
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  next3d(clear=F)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}

rgl.snapshot("pictures/MUScelltypes_long.png")
close3d()
}

#plot head and digestive muscles
{
plot_two_panel_background()
next3d(clear=F)
counter = 0
for (j in c(78:89)){    #iterate through the celltype list 
  counter = counter+1; print(j)
  #assign color
  color= color_head
  cells=celltypelist[[j]]
  #smooth with sigma 6000
  cells_smoothed=nlapply(cells, function(x) smooth_neuron(x, sigma=6000))
  #plot cells
  plot3d(cells_smoothed, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  next3d(clear=F)
  skids_raw=left_skids[j]
  #retrieve from catmaid all left cells of a type
  for (skids in skids_raw){
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
}

rgl.snapshot("pictures/MUScelltypes_head_dig.png")
close3d()
}

# make table of cells per segment per muscle type -------------------------

{
annotation_non_neuronal_celltypelist = list()
#read all non-neuronal celltypes - muscles from 37-89 and all annotations
for (i in c(37:89)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_non_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define the six body regions, matching the catmaid annotations
regions <- c('episphere','peristomium', 'segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')
sides <- c('left_side','right_side', 'middle')
types <- paste('celltype_non_neuronal', 37:89, sep='')

#we retrieve those skids that match our annotations
#three iterated lapply functions to retrieve skids by muscle type, body region and body side
muscle_per_side_per_body_region <- lapply(sides, function(s) lapply(regions, function(r) lapply(types, function(m) skids_by_annotation(annotation_non_neuronal_celltypelist,s,m,r))))

n_cells_left <- matrix(nrow=7,ncol=53)
n_cells_right <- matrix(nrow=7,ncol=53)
n_cells_middle <- matrix(nrow=7,ncol=53)
#count the occurrence of each type in each body region per side
for (j in 1:7){for (i in 1:53){
  n_cells_left[j,i] <- length(muscle_per_side_per_body_region[[1]][[j]][[i]])
  n_cells_right[j,i] <- length(muscle_per_side_per_body_region[[2]][[j]][[i]])
  n_cells_middle[j,i] <- length(muscle_per_side_per_body_region[[3]][[j]][[i]])
}}

library(heatmaply)
#add row and column names
row.names(n_cells_left) <- c('head left', 'peristomium left', 'sg0 left','sg1 left','sg2 left','sg3 left','pyg left')
row.names(n_cells_right) <- c('head right', 'peristomium right', 'sg0 right','sg1 right','sg2 right','sg3 right','pyg right')
row.names(n_cells_middle) <- c('head middle', 'peristomium middle', 'sg0 middle', 'sg1 middle','sg2 middle','sg3 middle','pyg middle')

#combine left and right matrix
n_cells <- rbind(n_cells_left, n_cells_right, n_cells_middle)

#celltype name lists
non_neuronal_celltype_names = c("akrotroch", "crescent cell", "prototroch", "nuchal cilia", "metatroch", "paratroch", "spinGland", "covercell", "ciliatedGland", "eyespot pigment cell", "pigment cell AE", "bright droplets parapodial", "macrophage", "yolk cover cell", "flat glia", "EC rad glia ", "MVGland", "microvillarCell", "protonephridium", "nephridium", "nephridiumTip", "chaeta", "aciculoblast", "circumacicular", "hemichaetal", "ER circumchaetal", "noER circumchaetal", "EC circumchaetal", "HeadGland", "InterparaGland", "spinMicroGland", "CB pigment", "vacuolar cell_head", "Glia pigmented", "pygidial pigment", "meso", "MUSac_notA", "MUSac_notP", "MUSac_notM", "MUSac_neuAV", "MUSac_neuPD", "MUSac_neuPV", "MUSac_neuDy", "MUSac_neuDx", "MUSac_neuDach", "MUSac_neure", "MUSac_i", "MUSob-ant_re", "MUSob-ant_arc", "MUSob-ant_m-pp", "MUSob-ant_ml-pp", "MUSob-ant_l-pp", "MUSob-ant_trans", "MUSob-post_notD", "MUSob-post_neuDlong", "MUSob-post_neuDprox", "MUSob-post_neuDdist", "MUSob-post_neuV", "MUSob-post_notV", "MUSob-postM", "MUSob-post_noty", "MUSob-post_i", "MUSchae_notDob", "MUSchae_notD", "MUSchae_notDn", "MUSchae_notA", "MUSchae_notAac", "MUSchae_notAre", "MUSchae_neuVob", "MUSchae_neuDac", "MUSchae_neuAVo", "MUSchae_neuAVt", "MUSchae_Are", "MUStrans_pyg", "MUSlong_D", "MUSlong_V", "MUSax", "MUSring", "MUSph", "MUSll", "MUSant", "MUSly", "MUSpl", "MUSci", "MUSch", "MUSpx", "MUSpr", "MUStri", "MUSmed_head", "EC")

#get the neuron name of the first skid in the skids list
non_neuronal_celltype_names <- list()
first_skids <- lapply(annotation_non_neuronal_celltypelist, function(x) x$skid[1])

for(i in 1:53){
  name<-catmaid_get_neuronnames(unlist(first_skids)[i], pid=11)
  name <- sub("_.*$", "", name )
  non_neuronal_celltype_names[i] <- name
}
length(non_neuronal_celltype_names)
#these are the muscle cell type names
colnames(n_cells)<-non_neuronal_celltype_names


#remove rows sg 0-3 middle (empty)
n_cells <- n_cells[-c(16:20),] 

as.data.frame(n_cells) %>%
  pivot_longer(cols = everything(), names_to = 'cell_type', values_to = 'number')

#convert to tibble
n_cells.tb <- rownames_to_column(as.data.frame(n_cells)) %>%
  pivot_longer(cols = starts_with("MUS"), 
               names_to = 'type', 
               values_to = 'number') %>%
  rename('region' = rowname)


  #geom_text(data=SN_IN_MN %>% filter(total_synapses>200), # Filter data first
   #         aes(label=Neuron), size=2.5, alpha=0.7, nudge_x =-0.1, check_overlap = TRUE, col='black')+
  
n_cells.tb %>%
  ggplot(aes(type, region, color = number)) +
  geom_point(shape = 15, stroke = 0, size = 4.7) +
  geom_text(data=n_cells.tb %>% filter(number>0), # Filter data first 
            aes(label = number), color = 'black', size = 3.2) +
  theme_half_open() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text (angle = 90, hjust = 1, vjust = 0.5, 
                                    size = 10, color = "black"), 
        axis.text.y = element_text (angle = 0, hjust = 1, vjust = 0.5, 
                                    size = 10, color = "black"),
        axis.title = element_blank(),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.7, 'mm'),
        axis.line = element_line(size = 0.3),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        ) +
  scale_colour_gradientn(
    colours = c('white', Okabe_Ito[c(1,6,7,5)]),
    na.value = "white",
    guide = "colourbar",
    guide_legend(title = '# of cells'),
    aesthetics = "colour") +
  scale_y_discrete(limits = c(rownames(n_cells)[1:7], rownames(n_cells)[16], rev(rownames(n_cells)[8:14]), rownames(n_cells)[15])) +
  #scale_y_discrete(limits = rev(rownames(n_cells))) +
  scale_x_discrete(limits = as.character(colnames(n_cells)))

# Saving plot
ggsave("pictures/Muscle_Celltypes_numbers.png", 
       width = 3000, height = 1150, limitsize = TRUE, 
       units = c("px"), bg = 'white')

#save table
write.csv(n_cells, file = "source_data/Figure2I_source_data_1.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


}

# tables of muscle names --------------------------------------------------

{
library(gtable)
  
not_ac <- tibble(subtype = c("", "", "Notopodial", 
                             "Neuropodial", "", "", "", "", "", "", 
                             "Interacicular"),
                 `cell type` = c("anterior notopodial acicular", 
                                       "posterior notopodial acicular", 
                                       "middle notopodial acicular",
                                       "anterior ventral neuropodial acicular",
                                       "posterior dorsal neuropodial acicular",
                                       "posterior ventral neuropodial acicular",
                                       "neuropodial Y",
                                       "dorsal neuropodial to notopodium",
                                       "dorsal neuropodial chaetal",
                                       "dorsal neuropodial acicular",
                                       "interacicular"))


not_ant_ob <- tibble(subtype = c("", "", "Ventral", 
                             "", "", ""),
                 `cell type` = c("parapodial retractor", 
                                 "ventral parapodial arc", 
                                 "medial oblique to mid-parapodium",
                                 "mediolateral oblique to mid-parapodium",
                                 "lateral oblique to mid-parapodium",
                                 "oblique to start of transverse"))

not_post_ob <- tibble(subtype = c("", "", "","Dorsal", 
                                 "Ventral", "", "", "", ""),
                     `cell type` = c("notopodial dorsal oblique", 
                                     "neuropodial dorsal oblique long", 
                                     "neuropodial dorsal oblique proximal",
                                     "neuropodial dorsal oblique distal",
                                     "posterior ventral neuropodial",
                                     "posterior ventral notopodial",
                                     "oblique to distal interacicular",
                                     "oblique to body wall",
                                     "distal interacicular"))

not_chae <- tibble(subtype = c("", "", "","", "","Notopodial", 
                                  "Neuropodial", "", "", "", ""),
                      `cell type` = c("notochaetal next to dorsal oblique", 
                                      "dorsal notopodial chaetal sac", 
                                      "next to dorsal notopodial chaetal sac",
                                      "anterior notopodial chaetal sac",
                                      "notopodial chaetal under acicula",
                                      "notopodial retractor",
                                      "neurochaetal next to ventral oblique",
                                      "neuropodial chaetal under acicula",
                                      "anterior ventral neurochaetal ob",
                                      "anterior ventral neurochaetal trans",
                                      "chaetal sac under parapodial retractor"))

not_head_dig <- tibble(subtype = c("", "Digestive", "", 
                               "", "Appendage", "", "",
                               "", "","other","", ""),
                   `cell type` = c("circular pygidial", 
                                   "pharyngeal", 
                                   "lower lip",
                                   "antenna",
                                   "lyrate",
                                   "palp",
                                   "cirrus",
                                   "cheek",
                                   "plexus",
                                   "ventral transverse of prostomium",
                                   "triangle", 
                                   "smooth head"
                                   ))

#define theme for tableGrob table
theme_1 <- ttheme_minimal(core = list(fg_params = list(hjust = 0, 
                                                       x = 0.03, 
                                                       fontsize = 5)),
                          colhead = list(fg_params = list(fontsize = 6, 
                                                          fontface = "bold")),
                          rowhead = list(fg_params = list(fontsize = 0)),
                          padding = unit(c(1, 1), "mm"))

#plot table and add grobs (horizontal lines)
table_ac <- tableGrob(not_ac, rows = NULL, theme = theme_1)
table_ac <- table_ac %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 4, b = 2, l = 1, r = ncol(table_ac)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 11, b = 2, l = 1, r = ncol(table_ac)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(1,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(1,"npc"),
    gp = gpar(lwd = 1)),
    t = 2, b = 2, l = 1, r = ncol(table_ac))
grid.newpage()
grid.draw(table_ac)

#plot table and add grobs (horizontal lines)
table_ant_ob <- tableGrob(not_ant_ob, rows = NULL, theme = theme_1)
table_ant_ob <- table_ant_ob %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(1,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(1,"npc"),
    gp = gpar(lwd = 1)),
    t = 2, b = 2, l = 1, r = ncol(table_ant_ob)) 
grid.newpage()
grid.draw(table_ant_ob)

#plot table and add grobs (horizontal lines)
table_post_ob <- tableGrob(not_post_ob, rows = NULL, theme = theme_1)
table_post_ob <- table_post_ob %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(1,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(1,"npc"),
    gp = gpar(lwd = 1)),
    t = 2, b = 2, l = 1, r = ncol(table_post_ob)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 5, b = 2, l = 1, r = ncol(table_post_ob)) 
grid.newpage()
grid.draw(table_post_ob)


#plot table and add grobs (horizontal lines)
table_chae <- tableGrob(not_chae, rows = NULL, theme = theme_1)
table_chae <- table_chae %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(1,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(1,"npc"),
    gp = gpar(lwd = 1)),
    t = 2, b = 2, l = 1, r = ncol(table_chae)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 7, b = 2, l = 1, r = ncol(table_chae)) 
grid.newpage()
grid.draw(table_chae)


#plot table and add grobs (horizontal lines)
table_head_dig <- tableGrob(not_head_dig, rows = NULL, theme = theme_1)
table_head_dig <- table_head_dig %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(1,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(1,"npc"),
    gp = gpar(lwd = 1)),
    t = 2, b = 2, l = 1, r = ncol(table_head_dig)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 4, b = 2, l = 1, r = ncol(table_head_dig)) %>%
  gtable_add_grob(grobs = segmentsGrob( # line across the bottom
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 1)),
    t = 8, b = 2, l = 1, r = ncol(table_head_dig)) 
grid.newpage()
grid.draw(table_head_dig)
}

# assemble figure ---------------------------------------------------------

{
panel_A <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_ac.png")) +
  draw_label("acicular", x = 0.5, y = 1.1, size = 9) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.43, y = 0.11, size = 6)  +
  geom_segment(aes(x = 0.45,
                   y = 0.93,
                   xend = 0.45,
                   yend = 0.82),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.6, "mm"))) +
  geom_segment(aes(x = 0.45,
                   y = 0.82,
                   xend = 0.45,
                   yend = 0.93),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.6, "mm"))) + 
  draw_label("a", x = 0.45, y = 0.96, size = 5) +
  draw_label("p", x = 0.45, y = 0.79, size = 5)  +
  geom_segment(aes(x = 0.55,
                   y = 0.93,
                   xend = 0.65,
                   yend = 0.93), 
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.6, "mm"))) +
  geom_segment(aes(x = 0.65,
                   y = 0.93,
                   xend = 0.55,
                   yend = 0.93),  
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.6, "mm"))) + 
  draw_label("v", x = 0.53, y = 0.93, size = 5) +
  draw_label("d", x = 0.67, y = 0.93, size = 5) 

panel_B <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_ant_ob.png")) +
  draw_label("anterior oblique", x = 0.5, y = 1.1, size = 9)
panel_C <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_post_ob.png")) +
  draw_label("posterior oblique", x = 0.5, y = 1.1, size = 9)
panel_D <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_chae.png")) +
  draw_label("chaetal", x = 0.5, y = 1.1, size = 9)

panel_E <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_trans.png")) +
  draw_label("transverse", x = 0.5, y = 1.1, size = 9)
panel_F <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_long.png")) +
  draw_label("longitudinal", x = 0.5, y = 1.1, size = 9) +
  draw_label("dorsolateral", x = 0.8, hjust = 0, y = 0.24, size = 5) +
  draw_label("ventrolateral", x = 0.8, hjust = 0, y = 0.19, size = 5) +
  draw_label("axochord", x = 0.8, hjust = 0, y = 0.14, size = 5)

panel_G <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_head_dig.png")) +
  draw_label("head and digestive tract", x = 0.5, y = 1.1, size = 9)
panel_H <- ggdraw() + draw_image(readPNG("pictures/MUScelltypes_all.png")) +
  draw_label("all muscles", x = 0.5, y = 0.98, size = 9)
panel_I <- ggdraw() + draw_image(readPNG("pictures/Muscle_Celltypes_numbers.png")) 

layout <- "
AAABB#CCCDD
EEEFF#GGGHH
IIIJJJKKKLL
MMMMNNNNNNN"

Figure2 <- panel_A + table_ac + panel_B + table_ant_ob + 
  panel_C + table_post_ob + panel_D + table_chae +
  panel_E + panel_F + panel_G + table_head_dig + 
  panel_H + panel_I +
  plot_layout(design = layout, heights = c(1, 1, 1, 1.6)) +
  plot_annotation(tag_levels = list(c('A', '', 'B', '',
                                      'C', '', 'D', '',
                                      'E', 'F', 'G', '', 
                                      'H', 'I'))) & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("figures/Figure2.pdf", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2600, height = 2400)


ggsave("figures/Figure2.png", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2600, height = 2400, bg='white')
}

