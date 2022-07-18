#Code to generate the anatomical panels of Figure 3 showing the different Leiden modules of the desmosomal connectome graph from Jasek et al 2021
#Code to generate morphological renderings for Figure3 of the Jasek et al. Desmosomal connectome paper
#needs as input a json file exported from gephi with module annotations (clusters)
#Gaspar Jekely 

source("code/Packages_and_Connection.R")

#load json file exported from gephi with the Leiden clusters colored
#would also work on a json file from catmaid but that uses HEX color codes so the relevant lines in the code would need to be skipped
#the Figure3-Leiden-modules.json file is also on github
graph2_json <- fromJSON(file = "data/Figure3-Leiden-modules.json")

cells_with_color <- data.table(id=numeric(), color=character())
setnames(cells_with_color, c("id"), c("skid"))

graph2_json$nodes[[1]]$attributes$unit_id
graph2_json$nodes[[1]]$color

cells_with_color[,1]
#change to data frame
setDF(cells_with_color)


#import from gephi json
for (i in c(1:length(graph2_json$nodes))){
  cells_with_color[i,1] <- graph2_json$nodes[[i]]$attributes$unit_id
  cells_with_color[i,2] <- graph2_json$nodes[[i]]$color
}
cells_with_color[1,]
#change back to data table
setDT(cells_with_color)
#set the key to be used for indexing
setkey(cells_with_color, color)
#to identify the column(s) indexed by
key(cells_with_color)
#order by color
setorder(cells_with_color, color)
skids <- as.numeric(unlist(cells_with_color[,1]))

#change rgb definiion by adding maxColorValue
cells_with_color[, color:=gsub(color,pattern=")",replacement = ",maxColorValue = 255)")]
cells_with_color[1:2,1:2]
#translate the rgb colors into HEX colors
#cells_with_color[,color:=eval(parse(text=color))]
cells_with_color[1:2,1:2]
colorHEX=list(1:nrow(cells_with_color))
for (i in c(0:nrow(cells_with_color))){print (eval(parse(text=cells_with_color[i,color])))
  colorHEX[i] <- eval(parse(text=cells_with_color[i,color]))}

#the partitions are named by their unique colors
partitions <- unique(cells_with_color[,2])
partitions
nrow(partitions)

# plotting by partition ---------------------------------------------------

#load anatomical references
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"), invertFaces = T, conn = NULL, pid = 11)
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"), invertFaces = T, conn = NULL, pid = 11)
acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11, fetch.annotations = F), function(x) smooth_neuron(x, sigma=6000))

## Function for desaturating colors by specified proportion
desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])}
#usage cc75 <- desat(cc, 0.75)

nopen3d() # opens a pannable 3d window

#iterate through partitions, read neurons from catmaid, plot and save png
for (i in c(1:nrow(partitions))){
  #retrieve a skid based on color key
  subset_skids <- cells_with_color[color==partitions[i],skid]
  
  mfrow3d(1, 1)  #defines the two scenes
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  par3d(zoom=0.52)
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.02,
         col="#E2E2E2",plotengine = getOption("nat.plotengine") )
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.05,
         col="#E2E2E2") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=0.3,
         col="black") 
  
  #read neurons based on skids and smooth them
  neurons <- nlapply(read.neurons.catmaid(subset_skids, pid=11), function(x) smooth_neuron(x, sigma=6000))

    #sample colors with different saturation 
  colors=list(1:length(neurons))
  for (j in c(1:length(neurons))){colors[j] <- desat(eval(parse(text=partitions[i])),sample(c(4:10/10),1))}

  plot3d(neurons, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
        rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
        col=as.character(colors)) 
  rgl.snapshot(paste("pictures/desmosomal_cluster", i, ".png",sep="") )
  clear3d()
}









