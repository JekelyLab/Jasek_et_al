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




# assembple figure --------------------------------------------------------

panelSEM <- ggdraw() + draw_image(readPNG("pictures/Platynereis_SEM_inverted2.png")) + 
  draw_label("Platynereis", x = 0.3, y = 0.99, size = 10, fontface = "italic") + 
  draw_label("larva", x = 0.55, y = 0.99, size = 10, fontface = "plain")

panel_mus <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus.png")) + 
  draw_label("all muscles", x = 0.3, y = 0.99, size = 10) +
  draw_label("aciculae", x = 0.9, y = 0.28, size = 8) +
  draw_line(
    x = c(0.58, 0.8, 0.79),
    y = c(0.14, 0.28, 0.42), size = 0.3) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.82),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.82,
                   xend = 0.1,
                   yend = 0.9),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.93, size = 8) +
  draw_label("p", x = 0.1, y = 0.79, size = 8) 


panel_des <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus_desmosomes.png")) + 
  draw_label("all desmosomes", x = 0.4, y = 0.99, size = 10)

panel_des_partners <- ggdraw() + 
  draw_image(readPNG("pictures/desmosomal_connectome_mus_partners.png")) + 
  draw_label("desmosomal partners", x = 0.45, y = 0.99, size = 10)

panel_des_conn <- ggdraw() + 
  draw_image(readPNG("pictures/full-graph2-darker.png")) + 
  draw_label("desmosomal connectome", x = 0.45, y = 0.99, size = 10) + 
  draw_label("#desmosomes", x = 0.86, y = 0.92, size = 8) + 
  draw_label("1", x = 0.84, y = 0.88, size = 8) + 
  draw_label("83", x = 0.84, y = 0.84, size = 8) +
  draw_line(x = c(0.88, 0.95), y = c(0.88, 0.88), size = 0.3, color = 'grey') +
  draw_line(x = c(0.88, 0.95), y = c(0.84, 0.84), size = 0.8, color = 'grey')

panel_desmo_EM <- ggdraw() + 
  draw_image(readPNG("pictures/desmo_EM_Fig1.png"))

layout <- "
AAAAAABBBBBBCCCCCCDDDDDD
EEEEEEEEEEEEEEEEFFFFFFFF
"

Figure1 <- panelSEM + panel_mus + panel_des + panel_des_partners +
  panel_desmo_EM + panel_des_conn +
  plot_layout(design = layout, guides = 'collect', heights = c(1,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("figures/Figure1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000)


ggsave("figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 2000, bg = 'white')

