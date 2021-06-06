#R code to generate figure panels for Figure 2 figure supplement 2 in Jasek et al. 2021 Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely March 2021

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)
#set working directory
setwd('/work_dir/')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

#retrieve cells and smooth them with sigma 6000
chaeta = nlapply(read.neurons.catmaid("^celltype_non_neuronal22$", pid=11, 
                                      fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

acicula = nlapply(read.neurons.catmaid("^celltype_non_neuronal23$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

circumacicular<-  nlapply(read.neurons.catmaid("^celltype_non_neuronal24$", pid=11, fetch.annotations = T), 
                          function(x) smooth_neuron(x, sigma=6000))	

hemichaetal<-  nlapply(read.neurons.catmaid("^celltype_non_neuronal25$", pid=11, fetch.annotations = T), 
                          function(x) smooth_neuron(x, sigma=6000))	

ER_circumchaetal <-  nlapply(read.neurons.catmaid("^celltype_non_neuronal26$", pid=11, fetch.annotations = T), 
                       function(x) smooth_neuron(x, sigma=6000))	

noER_circumchaetal<-  nlapply(read.neurons.catmaid("^celltype_non_neuronal27$", pid=11, fetch.annotations = T), 
                       function(x) smooth_neuron(x, sigma=6000))	

EC_circumchaetal<-  nlapply(read.neurons.catmaid("^celltype_non_neuronal28$", pid=11, fetch.annotations = T), 
                              function(x) smooth_neuron(x, sigma=6000))	


#load volumes based on catmaid id
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                              invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)



#############################################
#3d plotting of neurons
#############################################

#help to pick some colors
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(9, 'Blues')

# 3d plotting of cells
nopen3d(); mfrow3d(1, 2)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(20, 30, 1200, 800)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1))
par3d(zoom=0.52)


#plot meshes and background reference cells
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey60")
plot3d(circumacicular, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FD8D3C")

#move to next panel in rgl window
next3d(clear=F)
#define view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.52)

#plot 
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey60")
plot3d(circumacicular, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FD8D3C")

#make a snapshot to the working directory
rgl.snapshot("acicula_circumacicular.png")

##################################################
#plot chaetae with circumchaetal

#help to pick some colors
brewer.pal(9, 'Set3')

#move to next panel in rgl window
next3d(clear=F)
#clear it
clear3d()
#plot meshes and background reference cells
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey60")
plot3d(hemichaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FD8D3C")
plot3d(ER_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#80B1D3")
plot3d(noER_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#B3DE69")
plot3d(EC_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FCCDE5")

#move to next panel in rgl window
next3d(clear=F)
#define view
nview3d("right", extramat=rotationMatrix(pi, 0, 1, 1))
#define a sagittal clipping plane and re-zoom
clipplanes3d(1, 0, 0.16, -75700)
par3d(zoom=0.52)

plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey60")
plot3d(hemichaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FD8D3C")
plot3d(ER_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#80B1D3")
plot3d(noER_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#B3DE69")
plot3d(EC_circumchaetal, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="#FCCDE5")


#make a snapshot to the working directory
rgl.snapshot("chaeta_circumchaetal.png")


#######################################################
#display desmosomes

# 3d plotting of cells
nopen3d(); mfrow3d(1, 1)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(20, 30, 400, 198)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1))
par3d(zoom=0.16)
clipplanes3d(1, 0, 0.16, -75700)
rgl.bg(color='grey80')

#plot aciculae
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey10")
plot3d(circumacicular, WithConnectors = T, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#FD8D3C")
rgl.snapshot("acicula_circumacicular_desmosomes_closeup.png")

clear3d()

#plot aciculae as reference, but without connectors
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey10")

#plot chaetae and others with connectors
plot3d(chaeta, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col="grey60")
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
clipplanes3d(1, 0, 0.16, -79700)

plot3d(hemichaetal, WithConnectors = T, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#FD8D3C")
plot3d(ER_circumchaetal, WithConnectors = T, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#80B1D3")
plot3d(noER_circumchaetal, WithConnectors = T, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#B3DE69")
plot3d(EC_circumchaetal, WithConnectors = T, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="#FCCDE5")


rgl.snapshot("chaeta_circumchaetal_desmosomes_closeup.png")



###################################################
#we retrieve and plot the other cell types and basal lamina from the desmosomal connectome that connect to muscles

annot_desm_non_muscle <- catmaid_get_annotations_for_skeletons('^desmosome_connectome_non_muscle$', pid = 11)
head(annot_desm_non_muscle)
nlapply(read.neurons.catmaid("^bounding_dots$", pid=11, 
                             fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
cilia_with_desm <- nlapply(read.neurons.catmaid(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='ciliated cell'],
                                                pid=11), function(x) smooth_neuron(x, sigma=6000))
EC_with_desm <- nlapply(read.neurons.catmaid(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='celltype_non_neuronal90'], 
                                             pid=11), function(x) smooth_neuron(x, sigma=6000))
glia_with_desm <- nlapply(read.neurons.catmaid(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='glia'], 
                                               pid=11), function(x) smooth_neuron(x, sigma=6000))
bl_with_desm <- nlapply(read.neurons.catmaid(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='basal lamina'], 
                                             pid=11), function(x) smooth_neuron(x, sigma=6000))
pigment_with_desm <- nlapply(read.neurons.catmaid(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='pigment cell'], 
                                                  pid=11), function(x) smooth_neuron(x, sigma=6000))

length(annot_desm_non_muscle$skid[annot_desm_non_muscle$annotation=='ciliated cell'])

#plotting
# 3d plotting of cells
nopen3d(); mfrow3d(1, 1)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(20, 30, 600, 800)); nview3d("ventral", extramat=rotationMatrix(pi/20, 0, 1, 1))
par3d(zoom=0.52)
rgl.bg(color='white')
#plot ciliated cells from desmosomal connectome
{
  clear3d()
  #plot meshes and background reference cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
       col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
       col="#E2E2E2") 
  plot3d(cilia_with_desm, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(length(cilia_with_desm), 'Reds'))
  rgl.snapshot("cilia_with_desm.png")
}

#plot glia cells from desmosomal connectome
{
  clear3d()
  #plot meshes and background reference cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
         col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
         col="#E2E2E2") 
  plot3d(glia_with_desm, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
       col=hcl.colors(length(glia_with_desm), 'Blues'))
  rgl.snapshot("glia_with_desm.png")
}



#plot EC cells from desmosomal connectome
{
  clear3d()
  #plot meshes and background reference cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
         col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
         col="#E2E2E2") 
  plot3d(EC_with_desm, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(EC_with_desm), 'Purples 3'))
  rgl.snapshot("EC_with_desm.png")
}

#plot pigment cells from desmosomal connectome
{
  clear3d()
  #plot meshes and background reference cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
         col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
         col="#E2E2E2") 
  plot3d(pigment_with_desm, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(pigment_with_desm), 'Teal'))
  rgl.snapshot("pigment_with_desm.png")
}


#plot basal lamina from desmosomal connectome
{
  clear3d()
  #plot meshes and background reference cells
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
        rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
         col="#E2E2E2") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.07,
         col="#E2E2E2") 
  plot3d(bl_with_desm, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(bl_with_desm), 'Spectral'))
  rgl.snapshot("bl_with_desm.png")
}
