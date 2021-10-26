#Natverse code to display and work with synapse and desmosomal connectors separately
#Gaspar Jekely, Sept 2021

library(natverse)

# catmaid connection
source("~/R/conn.R")

#read in dorsal MUSlong cells with connectors
MUSlong = nlapply(read.neurons.catmaid("^celltype_non_neuronal75$", pid=11, 
                                       fetch.annotations = T), function(x) smooth_neuron(x, sigma=6000))
attributes(MUSlong)
class(MUSlong)
names(MUSlong)
attr(MUSlong, "df")
attr(MUSlong, "names")

#we get the connectors
MUSlong_conn <- connectors(MUSlong)

#separate synapse (prepost 1) and desmosomal (prepost 3) connectors
MUSlong_conn1 <- MUSlong_conn[MUSlong_conn$prepost %in% c("1"),]
MUSlong_conn3 <- MUSlong_conn[MUSlong_conn$prepost %in% c("3"),]
str(MUSlong_conn1)

#plot synapse and desmosomal connectors separately
plot3d(MUSlong_conn1$x, MUSlong_conn1$y, MUSlong_conn1$z, add = TRUE, col = 'red', size=6, alpha=0.5)
plot3d(MUSlong_conn3$x, MUSlong_conn3$y, MUSlong_conn3$z, add = TRUE, col = 'cyan', size=6, alpha=0.5)

MUSlong_conn3$x

