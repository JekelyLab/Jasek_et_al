#Code to generate Video3 of the desmosomal connectome paper Jasek et al
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

#read the saved visNetwork file from supplements/
conn_graph.visn <- readRDS("source_data/Figure3_source_data2.rds")


# plot neurons by colours matching network pic -----------------------------

#get names and colours and make df
names <- sapply(conn_graph.visn$nodes$id, function(x) x)
colors <- sapply(conn_graph.visn$nodes$color, function(x) x)
module <- sapply(conn_graph.visn$nodes$partition, function(x) x)
df <- data.frame(names, colors, module)
  
#get skid for module 4, 2nd cell
df[df$module==4,][2,][1]
  
#iterate through the neurons and plot them in the colour that was used in the network plot
plot_background_ventral()  

for (i in 1:13){
      print(i)
  filename = paste("videos/Video3_desmo_conn_module_", i+100, ".png", sep = "")
  rgl.snapshot(filename)
      for (j in 1:nrow(df[df$module==i,])){
        #read the skeleton based on the skid
        skeleton <- nlapply(read.neurons.catmaid(df[df$module==i,][j,][[1]], pid=11),
                            function(x) smooth_neuron(x, sigma=6000))
        #plot skeleton
        plot3d(skeleton, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
               rev = FALSE, fixup = F, add=T, alpha=1, col=df[df$module==i,][j,2])
        }

}

#full rotation
for (i in 100:427){
  play3d( spin3d( axis = c(0, 0, 10), rpm = 0.2), duration = 1)
  next3d(clear=F)
  print (i)
  #save a snapshot in the working directory
  rgl.snapshot(paste("videos/Video3_desmo_conn_module_spin_", i, ".png", sep = ""))
}
close3d()

#read png files and write video
library(av)
av::av_encode_video(paste('videos/', list.files("videos", '*.png'), sep = ""), 
                    framerate = 5,
                    output = 'videos/Video3.mp4')


file.remove(paste('videos/', list.files("videos", '*.png'), sep = ""))
#delete individual video frames




