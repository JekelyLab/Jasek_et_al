#R/Natverse code to export frames visualising the different muscle groups (outlines) in the left 2nd segment of NAOMI
#Jasek et al 2022 desmosomal connectome paper
#Gaspar Jekely Feb 2022

# load nat and all associated packages
library(natverse)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")


# catmaid connection, needs username, password AND token - weird!
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}

#a color-blind-friendly palette to use for skeleton coloring
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
pie(c(1,1,1,1,1,1,1,1), col=Okabe_Ito, labels = Okabe_Ito)


#read volumes and background structures
{
  outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                                invertFaces = T, conn = NULL, pid = 11)
  yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
  #these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
  bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  acicula_sg2l = nlapply(read.neurons.catmaid("^acicula_sg2l$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  chaeta_sg2l = nlapply(read.neurons.catmaid("^chaeta_sg2l$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  scalebar_50um_anterior = read.neurons.catmaid("^scalebar_50um_anterior$", pid=11)
}


#function to read several skeleton groups by up to five annotations and plot them
read_plot_neurons_2nd_seg <- function(annotation1, annotation2,annotation3,
                                      annotation4,annotation5){
  neuron1 = nlapply(read.neurons.catmaid(annotation1, pid=11),
                    function(x) smooth_neuron(x, sigma=1000))
  if (hasArg(annotation2)){
    neuron2 = nlapply(read.neurons.catmaid(annotation2, pid=11),
                      function(x) smooth_neuron(x, sigma=1000))
  }
  if (hasArg(annotation3)){
    neuron3 = nlapply(read.neurons.catmaid(annotation3, pid=11),
                      function(x) smooth_neuron(x, sigma=1000))
  }
  if (hasArg(annotation4)){
    neuron4 = nlapply(read.neurons.catmaid(annotation4, pid=11),
                      function(x) smooth_neuron(x, sigma=1000))
  }
  if (hasArg(annotation5)){
    neuron5 = nlapply(read.neurons.catmaid(annotation5, pid=11),
                      function(x) smooth_neuron(x, sigma=1000))
  }
  plot_background_2nd_seg()
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="#E69F00") 
  if (hasArg(annotation2)){
    plot3d(neuron2, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#CC79A7")
  }
  if (hasArg(annotation3)){
    plot3d(neuron3, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#CC79A7")
  }
  if (hasArg(annotation4)){
    plot3d(neuron4, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#009E73")
  }
  if (hasArg(annotation5)){
    plot3d(neuron5, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#D55E00") 
  }
}

#function to plot background and clip the plane to show sg2 left
plot_background_2nd_seg <- function(x){
  nopen3d() # opens a pannable 3d window
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.07,
         col="#E2E2E2") 
  plot3d(acicula_sg2l, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey50") 
  plot3d(chaeta_sg2l, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey80") 
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.52)
  nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 0))
  #z-axis clip
  clipplanes3d(0, 0, -1, 125000)
  #z-axis clip from top
  clipplanes3d(0, 0, 1, -50000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, -80700)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
}

#read and plot skeletons based on annotations (up to five with the read_plot_neurons_2nd_seg function)
read_plot_neurons_2nd_seg("^OUTLINE parapodial retractor muscle$","^OUTLINE anterior notopodial acicular muscle$", 
                          "^OUTLINE posterior notopodial acicular muscle$", "^OUTLINE ventral parapodial muscle arc$", 
                          "^OUTLINE oblique to distal interacicular$")


#add text labels in 3D
#adjust x,y,z coordinates
texts3d(100000,121000, 82000, text = "MUSobA-re", col='black', cex = 2)
texts3d(100000,69000, 85500, text = "MUSac-notA", col='black', cex = 2)
texts3d(100000,83000, 115000, text = "MUSac-notP", col='black', cex = 2)
texts3d(110000,120000, 100000, text = "MUSobA-arc", col='black', cex = 2)
texts3d(90000,115000, 115000, text = "MUSobP-M", col='black', cex = 2)
texts3d(77000,92000, 90000, text = "aciculae", col='black', cex = 2)
texts3d(87000,92000, 86000, text = "chaetae", col='black', cex = 2)

#export rotation by frame for video
for (i in 1:240){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 2)
  print (i)
  #save a snapshot
  filename <- paste("pictures/Video_sg2l_Mus_Outlines1_spin", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
}

#read the next set of skeletons

#define labels

#export rotation


