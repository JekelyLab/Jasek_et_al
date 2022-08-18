#packages to source and details of CATMAID connection

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# load nat and all associated packages, incl catmaid
library(nat)
library(catmaid)

options(nat.plotengine = 'rgl')
require("graphics")

library(tidyverse)
library(tidygraph)
library(rjson)
library(data.table)

library(gridExtra) #to render tables as grobs (grid graphical objects)
library(grid)

library(networkD3)
library(igraph)
library(leiden)
library(rgexf)
library(visNetwork)
library(plotly)
library(autoimage) #to rotate the coordinates of the graph plots

library(beepr) #to run beep() after a section finished

library(cowplot)
library(patchwork)
library(magick)
library(png)
library(pdftools)

library(colorspace)   ## hsv colorspace manipulations
library(RColorBrewer)

library(hash)

#create directory for R-generated pictures for figure panels (ignored by git)
dir.create("pictures")
dir.create("figures")
# catmaid connection, needs username, password AND token - weird!
conn <- source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
#for the public server:
#catmaid_login(server="https://catmaid.jekelylab.ex.ac.uk/", authname="AnonymousUser")

#cb friendly color palette
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                 "#CC79A7", "#000000")
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
                 '#CC6677', '#882255', '#AA4499', '#DDDDDD')
M_Winding_Col <-  c('#00753F','#1D79B7','#5D8C90','#D4E29E','#FF8734','#E55560',
'#F9EB4D','#C144BC','#FF9CFF','#8C7700','#77CDFC','#FFDAC7','#E0B1AD','#9467BD',
'#D88052','#A52A2A','grey')

display.brewer.all(colorblindFriendly = TRUE)
brewer12 <- brewer.pal(12, 'Paired')
Reds <- brewer.pal(9, 'Reds')
pie(rep(1,6),col=Reds[3:9], Reds[3:9])
pie(rep(1,8),col=Okabe_Ito[1:8], Okabe_Ito[1:8])

#read volumes
#These volumes are 3D structures in the animal's body and provide background 
#for reconstructed neurons or help create images
{
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                                invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
acicula <-   nlapply(read.neurons.catmaid("^acicula$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
scalebar_50um_anterior = read.neurons.catmaid("^scalebar_50um_anterior$", pid=11)
scalebar_50um_ventral = read.neurons.catmaid("^scalebar_50um_ventral$", pid=11)

}

#plotting function for ventral view with yolk and acicula
plot_background_ventral <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col="#E2E2E2") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey70")
  par3d(zoom=0.48)
}


#plotting function for ventral view with yolk and acicula
plot_background_ventral_no_acic <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col="#E2E2E2") 
  par3d(zoom=0.48)
}

#function to retrieve skids based on two annotations
skids_by_2annotations <- function(annotation1,annotation2){
  annotations_cells = list()
  annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11,conn=conn_http1)
  #we retrieve those skeletons that are also annotated with the second annotation
  return(unlist(lapply(annotations_cells,function(x) x[x$annotation==annotation2,1])))
}


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


