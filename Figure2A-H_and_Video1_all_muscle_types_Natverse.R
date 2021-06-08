#Desmosomal connectome paper R - Natverse code to plot all muscle cell types
#code for Video1 and Figure 1 panels A-H in the Jasek et al 2021 desmosomal connectome paper
#Gaspar Jekely - 4th Feb 2021

# load nat and associated packages
library(natverse)
library(magick)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")

# catmaid connection, needs username, password AND token - weird!
# can run this in a separate file using source function  source("~/R/conn.R")
conn = catmaid_login(server="https://catmaid.jekelylab.ex.ac.uk/", authname="AnonymousUser")
setwd("/work_directory/")

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

#check some values
celltypelist[[38]]$`1450345`
annotation_celltypelist[[38]]$annotation

left_skids = list(); counter=0
for (df1 in annotation_celltypelist){    #iterate through the celltype list 
  counter=counter+1
  #add cells with left_side annotation to left_skids list
  left_skids[[counter]] <- df1[df1$annotation == "left_side",1]
}



##########################################
# 3d plotting 
##########################################

#define randomised color list for muscle groups
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
texts3d(35000,0, 0, text = "ventral view", col="black", cex = 1.5)

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
texts3d(115000,30000, 8020, text = "left view", col="black", cex = 1.5)
next3d(clear=F)


#define the lighting
clear3d(type = "lights"); 
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")
#set background color
bg3d("gray100")
next3d(clear=F)
rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")


#add celltype text
next3d(clear=F); 
texts3d(95000,0, 6500, text = "all muscles", col="black", cex = 1.5)
rgl.snapshot("MUScelltypes_all_background.png"); next3d(clear=F)

#these are the annotation numbers for the different muscle celltypes in catmaid, same in our lists
#acicular 37:47 - Blues 3
#chaetal 63:73 - OrRd
#oblique anterior 48:53 - Purple-Blue
#oblique posterior 54:62 - Emrld
#transverse 74 - brown
#longitudinal 75-77 - Terrain 2
#digestive-head 78-89 - Warm

counter = 0
for (j in c(37:89)){    #iterate through the celltype list 
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
    if(skids_raw=="integer(0)") break; #some muscle celltypes are in the middle of the body so do not have righ or left members
    cells = nlapply(read.neurons.catmaid(skids, pid=11, 
                                         fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
    plot3d(cells, soma=T, lwd=3, col=color[counter], add=T, alpha=1, forceClipregion = TRUE); bg3d(col="white")
  }
  #make snapshot
 rgl.snapshot(paste("MUScelltypes_all", j, ".png", sep = ""))
}


#light moves across the 3d scene
for (i in c(12:-12)){
  clear3d(type = "lights")
  next3d(clear=F)
  rgl.light(i*5, 30, diffuse = "gray70")
  rgl.light(i*5, 30, specular = "gray5")
  rgl.light(-i*5, -30, specular = "gray5")
  next3d(clear=F)
  rgl.light(i*5, 30, diffuse = "gray70")
  rgl.light(i*5, 30, specular = "gray5")
  rgl.light(-i*5, -30, specular = "gray5")
  rgl.snapshot(paste("MUScelltypes_light", -i+13, ".png", sep = ""))
}


#full rotation
for (i in 1:185){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  next3d(clear=F)
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  next3d(clear=F)
  print (i)
  #save a snapshot in the working directory
  rgl.snapshot(paste("MUScelltypes_spin", i, ".png", sep = ""))
}


