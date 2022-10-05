#R/Natverse code for Video 2 of the Jasek et al 2022 desmosomal connectome paper
#Gaspar Jekely, Sanja Jasek 2022

source("code/Packages_and_Connection.R")

# colorblind-friendly color palette
Okabe_ito_noblack <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load landmarks
yolk <- catmaid_get_volume(4, rval = "mesh3d", invertFaces = T, pid = 11)
acicula_sg2l = nlapply(read.neurons.catmaid("^acicula_sg2l$", pid=11), function(x) smooth_neuron(x, sigma=6000))
chaeta_sg2l = nlapply(read.neurons.catmaid("^chaeta_sg2l$", pid=11), function(x) smooth_neuron(x, sigma=6000))


# read and plot muscles
read_and_plot <- function (namepattern, segment, color, alpha, alpha_text) {
  xsg2l <- paste(namepattern, "_", segment, sep = "") # only segment 2 on the left
  xskids <- catmaid_query_by_name(xsg2l, pid=11) # include skeletons and outlines
  if (length(xskids)>0) {
    # add text labels
    if (length(txt_pos) > 1) {
      texts3d(unlist(txt_pos[[namepattern]]), 
              text = namepattern,
              col='black',
              cex = 2, alpha = alpha_text)
    }
    xskids <- xskids$skid
    neurons = nlapply(read.neurons.catmaid(xskids, pid=11),
                    function(y) smooth_neuron(y, sigma=1000))
    
    plot3d(neurons, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=3,
          rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=alpha,
          col=color)
    
  }
  # cleanup
  rm(xskids)
  rm(neurons)
  rm(namepattern)
  rm(xsg2l)
}


# define which muscles to plot
MUS_groups <- list(c("MUSobP-neuV", "MUSobP-neuDlong"),
                   "MUSac-neuAV",
                   "MUSac-neure",
                   c("MUSac-neuDach", "MUSac-neuPV"),
                   "MUSac-notM",
                   "MUSac-notA",
                   "MUSac-notP",
                   c("MUSac-i", "MUSac-neuDx"))

dir.create("videos/Video2")


# first rotation show all muscle groups ------------------------------------

txt_pos <- NULL

names_non_neuronal <- read.csv("data/non_neuronal_celltypes_names.csv")
names_MUS <- list()
for (i in c(37:89)) {
  annotation_MUS <- paste("celltype_non_neuronal", i, sep="")
  name_MUS <- names_non_neuronal %>% 
    filter(CATMAID.annotation == annotation_MUS) %>% 
    select(Name) %>%
    unlist() %>%
    unname()
  print(name_MUS)
  names_MUS[[name_MUS]] <- annotation_MUS
}


cols_obant = hcl.colors(16, palette='Purple-Blue')
color_obant = sample(cols_obant[1:6])
cols_obpost = hcl.colors(16, palette='Green-Yellow')
color_obpost = sample(cols_obpost[1:9])
cols_ac = hcl.colors(20, palette='Blues 3')
color_ac = sample(cols_ac[1:11])
cols_ch = hcl.colors(16, palette='OrRd')
color_ch = sample(cols_ch[1:11])
cols_tr = hcl.colors(10, palette='YlOrBr')
color_tr = sample(cols_tr[1:1])
color_MUS <- c(color_ac, color_obant, color_obpost, color_ch, color_tr)

nopen3d()

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
par3d(zoom=0.6)
nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 1))
#z-axis clip
clipplanes3d(0, 0, -1, 145000)
#z-axis clip from top
clipplanes3d(0, 0, 1, -66000)
#y-axis clip
clipplanes3d(1, 0, 0.001, -60000)
#x-axis clip
clipplanes3d(0.1, 1, 0, -57000)
par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view

# all sg2 volumes, including MUStrans
for (i in (1:38)) {
  MUS <- names(names_MUS)[i]
  read_and_plot(MUS, "sg2l", color_MUS[i], 1)
}

# sg3 MUStrans
read_and_plot("MUStrans", "sg3l", color_tr, 0.8)

#export rotation by frame for video
for (l in 1:90){
  #  play3d( spin3d( axis = c(0, 0, 10), rpm = 4), duration =0.1 )
  nview3d("ventral", extramat=rotationMatrix(0.2, 1, 1, 1) %*%rotationMatrix(pi*l/90, 0, 0, 1))
  print (l)
  #save a snapshot
  filename <- paste("./videos/Video2/Video2_sg2l_Mus_Outlines_spin",
                    "1",
                    "_all", "_frame",
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
}

# plot individual muscle groups -------------------------------------------
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


col <- 1
rotation <- 1
for (MUS_group in MUS_groups) {
  print(paste("rotation", rotation))
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
  #texts3d(77000,92000, 90000, text = "aciculae", col='black', cex = 2)
  #texts3d(77000,92000, 86000, text = "chaetae", col='black', cex = 2)
  
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.6)
  #z-axis clip
  clipplanes3d(0, 0, -1, 145000)
  #z-axis clip from top
  clipplanes3d(0, 0, 1, -66000)
  #y-axis clip
  clipplanes3d(1, 0, 0.001, -60000)
  #x-axis clip
  clipplanes3d(0.1, 1, 0, -57000)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  um <- par3d()$userMatrix
  
  for (MUS in MUS_group) {
    read_and_plot(MUS, "sg2l", Okabe_ito_noblack[col], 1, 1)
    if (col == length(Okabe_ito_noblack)) {
      col = 1
    }
    else {
      col = col+1
    }
  }
  
 print(paste("color index", col))
 # all sg2 volumes, including MUStrans
 for (j in (1:38)) {
   MUS <- names(names_MUS)[j]
   read_and_plot(MUS, "sg2l", color_MUS[j], 0.05, 0)
 }
 
 # sg3 MUStrans
 read_and_plot("MUStrans", "sg3l", color_tr, 0.05, 0)
 
 #export rotation by frame for video
 for (l in 1:30){
   # play3d( spin3d( axis = c(0, 0, 10), rpm = 4), duration =0.1 )
    #rotate in a loop (with l e.g. 1:90 for a 180 turn)
    nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 0, 1))
    #print (l)
    #save a snapshot
    filename <- paste("./videos/Video2/Video2_sg2l_Mus_Outlines_spin",
                      rotation,
                      MUS_group[1], "_frame",
                      formatC(l, digits = 2, flag = "0"),
                      ".png", sep = "")
    #print(filename)
    rgl.snapshot(filename)
  }
 rotation = rotation + 1 
}

close3d()




#read png files and write video
av::av_encode_video(paste('videos/Video2/', list.files("videos/Video2/", '*.png'), sep = ""), 
                    framerate = 12,
                    output = 'videos/Video2.mp4')

unlink("videos/Video2/", recursive = T)

