#R code to generate Video6 (related to Figure8-figure-supplement1) of the Jasek et al. Platynereis desmosomal connectome paper
#Uses Natverse and accesses the data on catmaid
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

# Motoneuron plotting -------------------------------------------------------------

MN_names <- c("MNbow", "MNacicX","MNchae", "vMN1-2","MNwave","MNbiramous", "MNhose", "MNcrab", "MNladder", "MNspider type1", "MNspider type2","MNring")
annotations_MN <- c("celltype67", "celltype69","celltype165","vMN1-2", "celltype68", "celltype63", "celltype66", "celltype65", "celltype151","celltype61","celltype62","celltype64")
annotations_MN_mus <- c("MNbow_mus", "MNacicX_MUS","MNche_mus","vMN1-2_mus","MNwave_mus","MNbiramous_mus", "MNhose_mus","MNcrab_mus", "MNladder_mus","MNspider_type1_mus",
                        "MNspider_type2_mus","MNring_mus")
annotations_MN_mus_des <- c("MNbow_mus_des", "MNacicX_MUS_des","MNche_mus_des","vMN1-2_mus_des","MNwave_mus_des", "MNbiramous_mus_des", "MNhose_mus_des","MNcrab_mus_des", "MNladder_mus_des",
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
  texts3d(47000,0, 0, text = MN_names[i], col= "Black", cex = 1.5)
  rgl.snapshot(paste("videos/Video6_A", i, MN_names[i], "_1.png", sep=""))
  #plot muscles clustered by nblast clusters
  plot3d(MN_mus, col=sample(Reds, length(MN_mus), replace = TRUE), 
         soma=F, add=T, lwd=3, alpha=0.7)
  #add a text label
  texts3d(47000,0, 7000, text = paste(MN_names[i], " muscle targets"), col= "Red", cex = 1.5)
  rgl.snapshot(paste("videos/Video6_A", i, MN_names[i], "_2.png", sep=""))
  plot3d(MN_mus_des, soma=T, lwd=2,
         alpha=0.6,
         col=sample(Greens, length(MN_mus_des), replace = TRUE))
  #add a text label
  texts3d(47000,0, 14000, text = "desmosomal targets", col=hcl.colors(1, palette='YlGn'), cex = 1.5)
  par3d(zoom=0.57)
  rgl.snapshot(paste("videos/Video6_A", i, MN_names[i], "_3.png", sep=""))
  
}

#plot all
clear3d() 
plot3d(yolk, lwd=2,
       add=T,alpha=0.05,
       col="#E2E2E2") 


for (i in c(1:12)){
  #read skeletons from catmaid by annotation
  MN = nlapply(read.neurons.catmaid(annotations_MN[i], pid=11, 
                                    fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus = nlapply(read.neurons.catmaid(annotations_MN_mus[i], pid=11, 
                                        fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  MN_mus_des  = nlapply(read.neurons.catmaid(annotations_MN_mus_des[i], pid=11, 
                                             fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))
  plot3d(MN, soma=T, lwd=4, WithConnectors = FALSE,
         add=T, alpha=1,
         col=sample(Blues, length(MN), replace = TRUE)) 
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1, rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1, col="white") 
  rgl.snapshot(paste("videos/Video6_B", i, MN_names[i], "_all1.png", sep=""))
  #plot muscles clustered by nblast clusters
  plot3d(MN_mus, col=sample(Reds, length(MN_mus), replace = TRUE), 
         soma=F, add=T, lwd=3, alpha=1)
  rgl.snapshot(paste("videos/Video6_B", i, MN_names[i], "_all2.png", sep=""))
  plot3d(MN_mus_des, soma=T, lwd=2,
         alpha=1,
         col=sample(Greens, length(MN_mus_des), replace = TRUE))
  par3d(zoom=0.57)
  rgl.snapshot(paste("videos/Video6_B", i, MN_names[i], "_all3.png", sep=""))
  
}

plot3d(outline, add=T,alpha=0.03,
       col="#E2E2E2") 

#full rotation
for (i in 100:257){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.1), duration = 2)
  next3d(clear=F)
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.1), duration = 2)
  next3d(clear=F)
  print (i)
  #save a snapshot in the working directory
  rgl.snapshot(paste("videos/Video6_C_MNall_spin", i, ".png", sep = ""))
}
close3d()


# assemble video ---------------------------------------------------------

#read png files and write video
library(av)
av::av_encode_video(paste('videos/', list.files("videos", '*.png'), sep = ""), 
                    framerate = 5,
                    output = 'videos/Video6.mp4')


#delete individual video frames
file.remove(paste('videos/', list.files("videos", '*.png'), sep = ""))








