#R code to generate Figure8-figure-supplement1 of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# Motoneuron plotting -------------------------------------------------------------

MN_names <- c("MNacicX","MNchae","vMN1-2","MNbow", "MNwave","MNbiramous", "MNhose", "MNcrab", "MNladder", "MNspider type1", "MNspider type2","MNring")
annotations_MN <- c("celltype69","celltype165","vMN1-2","celltype67", "celltype68", "celltype63", "celltype66", "celltype65", "celltype151","celltype61","celltype62","celltype64")
annotations_MN_mus <- c("MNacicX_MUS","MNche_mus","vMN1-2_mus","MNbow_mus", "MNwave_mus","MNbiramous_mus", "MNhose_mus","MNcrab_mus", "MNladder_mus","MNspider_type1_mus",
                        "MNspider_type2_mus","MNring_mus")
annotations_MN_mus_des <- c("MNacicX_MUS_des","MNche_mus_des","vMN1-2_mus_des","MNbow_mus_des", "MNwave_mus_des", "MNbiramous_mus_des", "MNhose_mus_des","MNcrab_mus_des", "MNladder_mus_des",
                            "MNspider_type1_mus_des","MNspider_type2_mus_des","MNring_mus_des")
colors <- hcl.colors(40, palette = "Blues 3", rev=F)
colors <- colors[15:40]
Blues <- brewer.pal(9, 'Blues')[4:9]
Reds <- brewer.pal(9, 'Reds')[3:8]
Greens <- brewer.pal(9, 'Greens')[3:8]

plot_background_ventral_no_acic()
par3d(zoom=0.57)

for (i in c(1:12)){
  clear3d() 
  #read skeletons from catmaid by annotation
  MN = nlapply(read.neurons.catmaid(annotations_MN[i], pid=11, 
                                          fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus = nlapply(read.neurons.catmaid(annotations_MN_mus[i], pid=11, 
                                              fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus_des  = nlapply(read.neurons.catmaid(annotations_MN_mus_des[i], pid=11, 
                                                   fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  plot3d(yolk, lwd=2,
       add=T,alpha=0.05,
       col="#E2E2E2") 
  plot3d(outline, add=T,alpha=0.03,
         col="#E2E2E2") 
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
  plot3d(MN, soma=T, lwd=4,
       add=T, alpha=1,
       col=sample(Blues, length(MN), replace = TRUE)) 
  #add a text label
  #texts3d(47000,0, 0, text = MN_names[i], col= "Black", cex = 1.5)
  rgl.snapshot(paste("pictures/", MN_names[i], "_1.png", sep=""))
  #plot muscles clustered by nblast clusters
  plot3d(MN_mus, col=sample(Reds, length(MN_mus), replace = TRUE), 
         soma=F, add=T, lwd=3, alpha=0.7)
#add a text label
#texts3d(47000,0, 7000, text = paste(MN_names[i], " muscle targets"), col= "Red", cex = 1.5)
  rgl.snapshot(paste("pictures/", MN_names[i], "_2.png", sep=""))
  plot3d(MN_mus_des, soma=T, lwd=2,
       alpha=0.6,
       col=sample(Greens, length(MN_mus_des), replace = TRUE))
#add a text label
#texts3d(47000,0, 14000, text = "desmosomal targets", col=hcl.colors(1, palette='YlGn'), cex = 1.5)
  if(i == 6){plot3d(scalebar_50um_ventral, lwd=2, color = "black")}
  par3d(zoom=0.57)
  rgl.snapshot(paste("pictures/", MN_names[i], "_3.png", sep=""))
  
}
close3d()

# assemble figure ---------------------------------------------------------

panel_MNbi <- ggdraw() + draw_image(readPNG("pictures/MNbiramous_3.png")) +
  draw_label("MNbiramous", x = 0.5, y = 0.98, size = 8) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
            x = 0.84, y = 0.15, size = 8) +
  geom_segment(aes(x = 0.1,
                   y = 0.9,
                   xend = 0.1,
                   yend = 0.82),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
  geom_segment(aes(x = 0.1,
                   y = 0.82,
                   xend = 0.1,
                   yend = 0.9),
               size = 0.3,
               arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
  draw_label("a", x = 0.1, y = 0.93, size = 7) +
  draw_label("p", x = 0.1, y = 0.80, size = 7) 

panel_MNwave <- ggdraw() + draw_image(readPNG("pictures/MNwave_3.png")) +
  draw_label("MNwave", x = 0.5, y = 0.98, size = 8)

panel_spi1 <- ggdraw() + draw_image(readPNG("pictures/MNspider type1_3.png")) +
  draw_label("MNspider-ant", x = 0.5, y = 0.98, size = 8)

panel_spi2 <- ggdraw() + draw_image(readPNG("pictures/MNspider type2_3.png")) +
  draw_label("MNspider-post", x = 0.5, y = 0.98, size = 8)

panel_hose <- ggdraw() + draw_image(readPNG("pictures/MNhose_3.png")) +
  draw_label("MNhose", x = 0.5, y = 0.98, size = 8)

panel_ladder <- ggdraw() + draw_image(readPNG("pictures/MNladder_3.png")) +
  draw_label("MNladder", x = 0.5, y = 0.98, size = 8)

panel_chae <- ggdraw() + draw_image(readPNG("pictures/MNchae_3.png")) +
  draw_label("MNchae", x = 0.5, y = 0.98, size = 8)

panel_ring <- ggdraw() + draw_image(readPNG("pictures/MNring_3.png")) +
  draw_label("MNring", x = 0.5, y = 0.98, size = 8)

panel_acicX <- ggdraw() + draw_image(readPNG("pictures/MNacicX_3.png")) +
  draw_label("MNacicX", x = 0.5, y = 0.98, size = 8)

panel_vMN12 <- ggdraw() + draw_image(readPNG("pictures/vMN1-2_3.png")) +
  draw_label("vMN1, vMN2", x = 0.5, y = 0.98, size = 8)

layout <- "
ABCDE
FGHIJ
"

Figure8_fig_suppl1 <- panel_MNbi + panel_MNwave + panel_spi1 + panel_spi2 + panel_hose +
  panel_ladder + panel_chae + panel_ring + panel_acicX + panel_vMN12 +
  plot_layout(design = layout, heights = c(1,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("figures/Figure8_figure_supplemen1.pdf", limitsize = FALSE, 
       units = c("px"), Figure8_fig_suppl1, width = 3000, height = 1600)

ggsave("figures/Figure8_figure_supplemen1.png", limitsize = FALSE, 
       units = c("px"), Figure8_fig_suppl1, width = 3000, height = 1600, bg='white')



