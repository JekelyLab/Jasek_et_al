#R code to generate figure panels for Figure 1 figure supplement 2 in Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely

source("code/Packages_and_Connection.R")

# retrieve cells and smooth them with sigma 6000 --------------------------

{
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

outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                              invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)
}

# 3d plotting -------------------------------------------------------------

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
rgl.snapshot("pictures/acicula_circumacicular.png")

# plot chaetae with circumchaetal -----------------------------------------

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
rgl.snapshot("pictures/chaeta_circumchaetal.png")


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
rgl.snapshot("pictures/acicula_circumacicular_desmosomes_closeup.png")

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
rgl.snapshot("pictures/chaeta_circumchaetal_desmosomes_closeup.png")
close3d()
close3d()


# retrieve and plot the other cell types and basal lamina from  --------

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
rgl.snapshot("pictures/cilia_with_desm.png")


#plot glia cells from desmosomal connectome
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
rgl.snapshot("pictures/glia_with_desm.png")




#plot EC cells from desmosomal connectome
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
rgl.snapshot("pictures/EC_with_desm.png")


#plot pigment cells from desmosomal connectome
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
rgl.snapshot("pictures/pigment_with_desm.png")



#plot basal lamina from desmosomal connectome
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
rgl.snapshot("pictures/bl_with_desm.png")
close3d()

# assemble figure --------------------------------------------------------

panel_acic_circ <- ggdraw() + draw_image(readPNG("pictures/acicula_circumacicular.png")) +  
  draw_label("aciculae", x = 0.4, y = 0.99, size = 10, hjust = 0) +
  draw_label("circumacicular", x = 0.4, y = 0.94, size = 10, hjust = 0, color = '#FD8D3C') +
  geom_segment(aes(x = 0.05,
                 y = 0.95,
                 xend = 0.05,
                 yend = 0.87),
             arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.05,
                   y = 0.87,
                   xend = 0.05,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.05, y = 0.98, size = 8) +
  draw_label("p", x = 0.05, y = 0.84, size = 8) +
  geom_segment(aes(x = 0.85,
                   y = 0.95,
                   xend = 0.93,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.93,
                   y = 0.95,
                   xend = 0.85,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("v", x = 0.82, y = 0.95, size = 8) +
  draw_label("d", x = 0.96, y = 0.95, size = 8) 

panel_chae_circ <- ggdraw() + 
  draw_image(readPNG("pictures/chaeta_circumchaetal.png")) + 
  draw_label("chaetae", x = 0.4, y = 0.99, size = 10, hjust = 0) +
  draw_label("hemichaetal", x = 0.4, y = 0.95, size = 10, 
             hjust = 0, color = '#FD8D3C') +
  draw_label("ER circumchaetal", x = 0.4, y = 0.91, size = 10, 
             hjust = 0, color = '#80B1D3') +
  draw_label("noER circumchaetal", x = 0.4, y = 0.87, size = 10, 
             hjust = 0, color = '#B3DE69') +
  draw_label("EC circumchaetal", x = 0.4, y = 0.83, size = 10, 
             hjust = 0, color = '#FCCDE5')  +
  draw_label("neuropodia", x = 0.54, y = 0.68, size = 10)  +
  draw_label("notopodia", x = 0.88, y = 0.64, size = 10)  +
  geom_segment(aes(x = 0.05,
                   y = 0.95,
                   xend = 0.05,
                   yend = 0.87),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.05,
                   y = 0.87,
                   xend = 0.05,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.05, y = 0.98, size = 8) +
  draw_label("p", x = 0.05, y = 0.84, size = 8) +
  geom_segment(aes(x = 0.85,
                   y = 0.95,
                   xend = 0.93,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.93,
                   y = 0.95,
                   xend = 0.85,
                   yend = 0.95),
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("v", x = 0.82, y = 0.95, size = 8) +
  draw_label("d", x = 0.96, y = 0.95, size = 8) +
  draw_line(x = c(0.53, 0.63), y = c(0.65, 0.62), size = 0.2) +
  draw_line(x = c(0.53, 0.62), y = c(0.65, 0.5), size = 0.2) +
  draw_line(x = c(0.53, 0.61), y = c(0.65, 0.3), size = 0.2) +
  draw_line(x = c(0.87, 0.76), y = c(0.6, 0.62), size = 0.2) +
  draw_line(x = c(0.87, 0.75), y = c(0.6, 0.5), size = 0.2) +
  draw_line(x = c(0.87, 0.74), y = c(0.6, 0.3), size = 0.2) +
  draw_line(x = c(0.2, 0.45, 0.45, 0.2, 0.2), 
            y = c(0.57, 0.57, 0.41, 0.41, 0.57), size = 0.3, linetype = 2)
  


panel_acic_circ_close <- ggdraw() + 
  draw_image(readPNG("pictures/acicula_circumacicular_desmosomes_closeup.png")) + 
  draw_label("aciculae", x = 0.97, y = 0.8, hjust = 1, size = 10) +
  draw_label("circumacicular", x = 0.97, y = 0.72, size = 10, 
       hjust = 1, color = '#FD8D3C')  +
  draw_label("desmosomes", x = 0.97, y = 0.64, size = 10, 
             hjust = 1, color = '#FFFFFF') 
panel_chae_circ_close <- ggdraw() + 
  draw_image(readPNG("pictures/chaeta_circumchaetal_desmosomes_closeup.png")) + 
  draw_label("aciculae", x = 0.97, y = 0.82, hjust = 1, size = 10) +
  draw_label("chaetae", x = 0.97, y = 0.74, size = 10, 
             hjust = 1, color = 'grey40')  +
  draw_label("desmosomes", x = 0.97, y = 0.66, size = 10, 
             hjust = 1, color = '#FFFFFF') 

panel_close <- ggdraw(panel_acic_circ_close / panel_chae_circ_close)

panel_cilia <- ggdraw() + 
  draw_image(readPNG("pictures/cilia_with_desm.png")) + 
  draw_label("ciliary band cells", x = 0.45, y = 0.99, size = 10) 

panel_glia <- ggdraw() + 
  draw_image(readPNG("pictures/glia_with_desm.png")) + 
  draw_label("glia cells", x = 0.45, y = 0.99, size = 10) 

panel_EC <- ggdraw() + 
  draw_image(readPNG("pictures/EC_with_desm.png")) + 
  draw_label("epidermal cells", x = 0.45, y = 0.99, size = 10)
panel_pigment <- ggdraw() + 
  draw_image(readPNG("pictures/pigment_with_desm.png")) + 
  draw_label("pigment cells", x = 0.45, y = 0.99, size = 10)
panel_bl <- ggdraw() + 
  draw_image(readPNG("pictures/bl_with_desm.png")) + 
  draw_label("basal lamina fragments", x = 0.45, y = 0.99, size = 10)


layout <- "
AAAAAAABBBBBBBCCCCCC
AAAAAAABBBBBBBCCCCCC
EEEEFFFFGGGGHHHHIIII
"

Figure1_figSuppl2 <- panel_acic_circ + panel_chae_circ + panel_close +
  panel_cilia +  panel_EC + panel_glia + panel_pigment + panel_bl +
  plot_layout(design = layout, 
              guides = 'collect', 
              tag_level = 'new',
              heights = c(0.5, 0.5, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("figures/Figure1_figSuppl2.pdf", limitsize = FALSE, 
       units = c("px"), Figure1_figSuppl2, width = 3600, height = 2100)


ggsave("figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3600, height = 2100, bg = 'white')

