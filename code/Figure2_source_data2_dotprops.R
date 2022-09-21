#Desmosomal connectome paper R - code to generate dotprops source data
#Gaspar Jekely

source("code/Packages_and_Connection.R")

MUS = nlapply(
  read.neurons.catmaid(
    "^muscle$", 
    pid=11
    ), function(x) smooth_neuron(
      x, sigma=6000)
)


# convert neuron representations to "dotprops" ie point-vector format for
# nblast - /1e3 scales to microns instead of nm - required for nblast v2
MUS_dots=dotprops(MUS/1e3, k=5, resample=1)

plot3d(MUS_dots)

saveRDS(MUS_dots, file = "source_data/Figure2_source_data2.RDS")

