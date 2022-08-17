#R/Natverse code to export frames visualising the different muscle groups (outlines) in the left 2nd segment of NAOMI
#Jasek et al 2022 desmosomal connectome paper
#Gaspar Jekely Feb 2022

library(catmaid)
library(nat)
library(hash)
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


# read and plot muscles
read_and_plot <- function (x, color) {
  xsg2l <- paste(x, "_sg2l", sep = "") # only segment 2 on the left
  xskids <- catmaid_query_by_name(xsg2l, pid=11) # include skeletons and outlines
  xskids <- xskids$skid
  neurons = nlapply(read.neurons.catmaid(xskids, pid=11),
                  function(y) smooth_neuron(y, sigma=1000))
  plot3d(neurons, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
        rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
        col=color)
  # add text labels
  texts3d(unlist(txt_pos[[x]]), text = x, col='black', cex = 2)
  # cleanup
  rm(xskids)
  rm(neurons)
  rm(x)
  rm(xsg2l)
}


# x,y,z coordinates for muscle text labels
txt_pos <- hash(
  "MUSobA-re" = c(105000,110000, 82000),
  "MUSobA-arc" = c(110000,120000, 100000),
  "MUSobA-mpp" = c(94000,105000, 77000),
  "MUSobA-lpp" = c(90000,115000, 80000),
  "MUSobA-trans" = c(11000,110000, 80000),
  "MUSobP-neuD" = c(110000,78000, 107000),
  "MUSobP-neuDlong" = c(98000,75000, 116000),
  "MUSobP-neuV" = c(98000,120000, 111000),
  "MUSac-neuAV" = c(88000,114000, 90000),
  "MUSac-neure" = c(76000,100000, 87000),
  "MUSchae-neuDac" = c(104000,99000, 95000),
  "MUSac-neuDach" = c(90000,94000, 86000),
  "MUSac-neuPV" = c(86000,123000, 107000),
  "MUSac-neuDy" = c(92000,100000, 108000),
  "MUSac-notM" = c(90000,95000, 85000),
  "MUSac-notA" = c(92000,70000, 88000),
  "MUSchae-notAac" = c(100000,85000, 115000),
  "MUSac-notP" = c(100000,73000, 112000),
  "MUSac-i" = c(69000,95000, 87000),
  "MUSac-neuDx" = c(85000,92000, 103000)
)


# define which muscles to plot
MUS_groups <- list(c("MUSobP-neuV", "MUSobP-neuDlong"), "MUSac-neuAV", "MUSac-neure", c("MUSac-neuDach", "MUSac-neuPV"), "MUSac-notM", "MUSac-notA", "MUSac-notP", c("MUSac-i", "MUSac-neuDx"))

dir.create("pictures/Video-PPC")

nopen3d() # opens a pannable 3d window

i = 1
for (MUS_group in MUS_groups)
{
  clear3d()
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
  
  # add text labels for chaetae and aciculae
  texts3d(77000,92000, 90000, text = "aciculae", col='black', cex = 2)
  texts3d(77000,92000, 86000, text = "chaetae", col='black', cex = 2)
  
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.6)
  nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 1))
  #z-axis clip
  clipplanes3d(0, 0, -1, 133000)
  #z-axis clip from top
  clipplanes3d(0, 0, 1, -66000)
  #y-axis clip
  clipplanes3d(1, 0, 0.001, -60000)
  #x-axis clip
  clipplanes3d(0.1, 1, 0, -57000)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  for (MUS in MUS_group) {
    read_and_plot(MUS, Okabe_Ito[i])

    if (i == length(Okabe_Ito)) {
      i = 1
    }
    else {
      i = i+1
    }
  }

  #export rotation by frame for video
  for (l in 1:124){
    play3d( spin3d( axis = c(0, 0, 10), rpm = 2), duration =0.25 )
    print (l)
    #save a snapshot
    filename <- paste("./pictures/Video-PPC/Video_sg2l_Mus_Outlines_spin",  formatC(i, digits = 1, flag = "0"), MUS_group[1], "_frame", formatC(l, digits = 1, flag = "0"), ".png", sep = "")
    #print(filename)
    rgl.snapshot(filename)
  }

}
close3d()

#read png files and write video
library(av)
av::av_encode_video(paste('pictures/Video-PPC/', list.files("pictures/Video-PPC/", '*.png'), sep = ""), 
                    framerate = 5,
                    output = 'videos/Video-PPC.mp4')
