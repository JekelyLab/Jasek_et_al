#R code to generate panels B, C and D of Figure1 of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# read cells and connectors --------------------------------------------------------------
{
catmaid_get_volumelist(conn = NULL, pid = 11)
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11)

acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, 
                                      fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
muscle = nlapply(read.neurons.catmaid("^muscle$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))

desmosome_connectome_non_muscle = nlapply(read.neurons.catmaid("^desmosome_connectome_non_muscle$", pid=11, 
                                                                   fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))


#get the connectors
MUSconn <- connectors(muscle)

#separate synapse (prepost 1) and desmosomal (prepost 3) connectors
MUSconn1 <- MUSconn[MUSconn$prepost %in% c("1"),]
MUSconn3 <- MUSconn[MUSconn$prepost %in% c("3"),]

}

# 3d plotting -------------------------------------------------------------

{
#do nblast of muscles so that they are coloured by clusters
# nblast muscles - /1e3 scales to microns instead of nm - required for nblast v2
np_dots=dotprops(muscle/1e3, k=5, resample=1)
nblast_v2 = nblast_allbyall(np_dots)
# clustering nblast results: v2
v2_clust = nhclust(scoremat=nblast_v2, method="ward.D2") # method arg refers to
# hierarchical clustering method - default is ward.D, an incorrect version of
# Ward's algorithm - i prefer D2, or average linkage


nopen3d() # opens a pannable 3d window
mfrow3d(1, 1)  #defines the two scenes
par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
par3d(zoom=0.52)

#panel B of Figure 2
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 

plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 

plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 

plot3d(v2_clust, k=24, db=muscle, lwd=1,col=hcl.colors(45, palette='Reds'))

#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus.png")


#panel C of Figure 1
clear3d()
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 
plot3d(MUSconn3$x, MUSconn3$y, MUSconn3$z, add = TRUE, 
       col=hcl.colors(20000, palette='Oranges'), size=4, alpha=0.5)

#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus_desmosomes.png")

#panel D of Figure 2
clear3d()
plot3d(desmosome_connectome_non_muscle, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col=hcl.colors(2000, palette='Blues')) 
plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.04,
       col="#E2E2E2") 
plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.1,
       col="#E2E2E2") 
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.5,
       col="black") 
#make snapshot
rgl.snapshot("pictures/desmosomal_connectome_mus_partners.png")
close3d()
}

