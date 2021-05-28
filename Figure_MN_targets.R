#Desmosomal connectome paper R - Natverse code to plot all muscle targets of MNs and the desmosomal partners of those muscles
#code for Figures showing MNs and their Muscle and desmosomal targets Jasek et al 2021 paper
#Gaspar Jekely - 5th Feb 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# load nat and all associated packages, incl catmaid
library(natverse)
library(magick)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")
# catmaid connection, needs username, password AND token
# load nat and all associated packages, incl catmaid

# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
source("~/R/conn.R")
workdir <- "/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Muscles/Figures/Figure8-figure-supplement1-MN-targets"
#workdir <- "/your_working_directory/"
setwd(workdir)

catmaid_get_volumelist(conn = NULL, pid = 11)
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

chaeta = nlapply(read.neurons.catmaid("^chaeta$", pid=11, 
                                          fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, 
                                      fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11, 
                                                                   fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))


##########################################
#MN plotting
##########################################

MN_names <- c("MNsmile","MNche","vMN1-2","MNbow", "MNwave","MNbiramous", "MNhose", "MNcrab", "MNladder", "MNspider type1", "MNspider type2","MNring")
annotations_MN <- c("celltype182","celltype165","vMN1-2","celltype67", "celltype68", "celltype63", "celltype66", "celltype65", "celltype151","celltype61","celltype62","celltype64")
annotations_MN_mus <- c("MNsmile_mus","MNche_mus","vMN1-2_mus","MNbow_mus", "MNwave_mus","MNbiramous_mus", "MNhose_mus","MNcrab_mus", "MNladder_mus","MNspider_type1_mus",
                        "MNspider_type2_mus","MNring_mus")
annotations_MN_mus_des <- c("MNsmile_mus_des","MNche_mus_des","vMN1-2_mus_des","MNbow_mus_des", "MNwave_mus_des", "MNbiramous_mus_des", "MNhose_mus_des","MNcrab_mus_des", "MNladder_mus_des",
                            "MNspider_type1_mus_des","MNspider_type2_mus_des","MNring_mus_des")
colors=hcl.colors(40, palette = "Blues 3", rev=F)
length(MN_names)
nopen3d(); mfrow3d(1, 1)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 600, 800)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)
#define the lighting
clear3d(type = "lights"); rgl.light(60, 30, diffuse = "gray70"); rgl.light(60, 30, specular = "gray5"); rgl.light(-60, -30, specular = "gray5")


for (i in c(1:12)){
#read skeletons from catmaid by annotation
MN = nlapply(read.neurons.catmaid(annotations_MN[i], pid=11, 
                                          fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
MN_mus = nlapply(read.neurons.catmaid(annotations_MN_mus[i], pid=11, 
                                              fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
MN_mus_des  = nlapply(read.neurons.catmaid(annotations_MN_mus_des[i], pid=11, 
                                                   fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
# NBlast of muscles to color them by cluster
#first convert to dotprops for Nblast
MN_MUS_dots=dotprops(MN_mus/1e3, k=7, resample=1)
#nblast 
blast_scores <- nblast_allbyall(MN_MUS_dots)
# Hierarchically cluster the n.blast scores
hckcs <- nhclust(scoremat=blast_scores)
# 3d plotting of clustering results and other cells
#plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
#       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
#       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.15,
       col="#E2E2E2") 
plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
#plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=1, rev = FALSE, fixup = F, add=T, 
#       forceClipregion = TRUE, alpha=0.6, col="gray10") 
plot3d(MN, WithConnectors = F, NeuronNames = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=colors[1:length(MN)])
#add a text label
#texts3d(47000,0, 0, text = MN_names[i], col= "Black", cex = 1.5)
rgl.snapshot(paste(MN_names[i], "_1.png", sep=""))
#plot muscles clustered by nblast clusters
plot3d(hckcs, k=10, col=hcl.colors(15, palette='OrYel'), db=MN_mus, soma=F, add=T, lwd=3, alpha=1)
#add a text label
#texts3d(47000,0, 7000, text = paste(MN_names[i], " muscle targets"), col= "Red", cex = 1.5)
rgl.snapshot(paste(MN_names[i], "_2.png", sep=""))
plot3d(MN_mus_des, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.6,
       col=hcl.colors(length(MN_mus_des), palette='YlGn'))
#add a text label
#texts3d(47000,0, 14000, text = "desmosomal targets", col=hcl.colors(1, palette='YlGn'), cex = 1.5)
rgl.snapshot(paste(MN_names[i], "_3.png", sep=""))
clear3d() 
}




#plot all MNs and their targets without clearing the RGL window
for (i in c(1:12)){
  #read skeletons from catmaid by annotation
  MN = nlapply(read.neurons.catmaid(annotations_MN[i], pid=11, 
                                    fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus = nlapply(read.neurons.catmaid(annotations_MN_mus[i], pid=11, 
                                        fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus_des  = nlapply(read.neurons.catmaid(annotations_MN_mus_des[i], pid=11, 
                                             fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  # NBlast of muscles to color them by cluster
  #first convert to dotprops for Nblast
  MN_MUS_dots=dotprops(MN_mus/1e3, k=7, resample=1)
  #nblast 
  blast_scores <- nblast_allbyall(MN_MUS_dots)
  # Hierarchically cluster the n.blast scores
  hckcs <- nhclust(scoremat=blast_scores)
  # 3d plotting of clustering results and other cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
         col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
         col="#E2E2E2") 
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=1, rev = FALSE, fixup = F, add=T, 
         forceClipregion = TRUE, alpha=0.4, col="gray10") 
  plot3d(MN, WithConnectors = F, NeuronNames = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=colors[1:length(MN)])
  rgl.snapshot(paste(MN_names[i], "_all1.png", sep=""))
  #plot muscles clustered by nblast clusters
  plot3d(hckcs, k=10, col=hcl.colors(15, palette='OrYel'), db=MN_mus, soma=F, add=T, lwd=3, alpha=1)
    rgl.snapshot(paste(MN_names[i], "_all2.png", sep=""))
  plot3d(MN_mus_des, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.6,
         col=hcl.colors(length(MN_mus_des), palette='YlGn'))
  rgl.snapshot(paste(MN_names[i], "_all3.png", sep=""))
}



#light moves across the 3d scene
for (i in c(12:-12)){
  clear3d(type = "lights")
  rgl.light(i*5, 30, diffuse = "gray70")
  rgl.light(i*5, 30, specular = "gray5")
  rgl.light(-i*5, -30, specular = "gray5")
  rgl.snapshot(paste("MN_targets_light", -i+13, ".png", sep = ""))
}

nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1)); par3d(zoom=0.5)
#full rotation
for (i in 1:185){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  print (i)
  Sys.sleep(10)
  #save a snapshot in the working directory
  rgl.snapshot(paste("MN_targets_spin", i, ".png", sep = ""))
}


