#packages to source and details of CATMAID connection

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# load nat and all associated packages, incl catmaid
library(natverse)
library(magick)
library(nat)
options(nat.plotengine = 'rgl')
require("graphics")

library(networkD3)
library(igraph)
library(leiden)
library(beepr) #to run beep() after a section finished
library(cowplot)
library(pdftools)
library(tidyverse)
library(magick)
library(rjson)
library(data.table)
library(colorspace)   ## hsv colorspace manipulations

library(png)
library(patchwork)


#create directory for R-generated pictures for figure panels (ignored by git)
dir.create("pictures")
dir.create("figures")
# catmaid connection, needs username, password AND token - weird!
conn <- source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
#for the public server:
#catmaid_login(server="https://catmaid.jekelylab.ex.ac.uk/", authname="AnonymousUser")

#cb friendly color palette
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                 "#CC79A7", "#000000")
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
                 '#CC6677', '#882255', '#AA4499', '#DDDDDD')
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
brewer12 <- brewer.pal(12, 'Paired')



