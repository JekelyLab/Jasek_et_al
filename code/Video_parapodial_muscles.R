#R/Natverse code to export frames visualising the different muscle groups (outlines) in the left 2nd segment of NAOMI
#Jasek et al 2022 desmosomal connectome paper
#Gaspar Jekely Feb 2022

library(natverse)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")

# colorblind-friendly color palette
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# sourcing catmaid login, see https://natverse.org/rcatmaid/reference/catmaid_login.html
conn <- source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
#conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))

# load landmarks
yolk <- catmaid_get_volume(4, rval = "mesh3d", invertFaces = T, pid = 11)
acicula_sg2l = nlapply(read.neurons.catmaid("^acicula_sg2l$", pid=11), function(x) smooth_neuron(x, sigma=6000))
chaeta_sg2l = nlapply(read.neurons.catmaid("^chaeta_sg2l$", pid=11), function(x) smooth_neuron(x, sigma=6000))


nopen3d() # opens a pannable 3d window

# plot landmarks
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


# read and plot muscles

MUS_list <- c("MUSobA-re_sg2l","MUSac-notA_sg2l", "MUSac-notP_sg2l", "MUSobA-arc_sg2l", "MUSobP-M_sg2l")

read_and_plot <- function (x, color) {
  xskids <- catmaid_query_by_name(x, pid=11)
  xskids <- xskids$skid

  neurons = nlapply(read.neurons.catmaid(xskids, pid=11),
                  function(x) smooth_neuron(x, sigma=1000))
  plot3d(neurons, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
        rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
        col=color)
}


i = 1
while (i <= length(Okabe_Ito) && i <= length(MUS_list)) {
  read_and_plot(MUS_list[i], Okabe_Ito[i])
  i = i+1
}

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