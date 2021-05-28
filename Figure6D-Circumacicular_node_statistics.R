#weighted degrees (WD) and degrees (D) of nodes in the desmosomal connectome
#plot for Figure 5 of Jasek et al 2021 Desmosomal connectome paper
setwd("/Users/gaspar/OneDrive\ -\ University\ of\ Exeter/Paper/Muscles/Figures/Figure5-acicula-weights/")


WD_circumacicular <- c(83,51,50,50,48,47,45,43,42,41,40,40,39,39,38,38,36,36,34,34,33,33,31,31,31,30,29,29,29,28,25,23,22,21,21,20,19,19,18,17,17,15,15,15,14,12,12,12,12,10,8,4,3,2)
WD_MUSchae_notAac <- c(67,60,56,49,40,38,31,31,30,29,27,24,24,23,19,17,17,16,16,15,10,10,9,9)
WD_MUSchae_neuDac <- c(48,41,38,37,36,34,33,32,29,26,26,22,21,21,20,19,18,18,16,14,14,12,9,8,8,6,5,4)
WD_MUSac_notM_P_A <- c(57,56,40,35,34,34,33,31,29,28,28,27,27,26,25,22,21,21,20,20,20,18,17,17,16,12,11,9,9,8,8,7,7,6,6,6,4)
WD_acicula <- c(39, 39, 28, 18, 15, 14, 13, 12, 11, 10, 9, 7)
WD_MUSlong_V <- c(47, 42, 30, 29, 27, 26, 24, 24, 21, 21, 21, 21, 19, 19, 18, 18, 17, 16, 16, 16, 16, 16, 15, 15, 15, 14, 14, 14, 14, 13, 13, 13, 12, 12, 12, 11, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 2, 2, 2)
WD_MUSlong_D <- c(33, 27, 26, 25, 25, 24, 24, 23, 22, 22, 22, 21, 20, 20, 20, 20, 19, 18, 17, 16, 15, 15, 15, 14, 13, 13, 13, 13, 13, 13, 12, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2)
WD_circumchaetal <- c(20, 17, 17, 17, 16, 16, 15, 15, 14, 13, 13, 13, 13, 12, 12, 12, 12, 11, 11, 11, 10, 10, 10, 10, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2, 2, 1)

D_circumacicular <- c(9, 6, 8, 7, 10, 9, 7, 10, 8, 8, 9, 7, 13, 9, 8, 6, 6, 9, 6, 8, 7, 6, 7, 8, 7, 5, 5, 6, 5, 6, 12, 7, 7, 5, 8, 5, 6, 7, 4, 3, 5, 4, 3, 5, 4, 6, 4, 4, 2, 4, 4, 3, 2, 2)
D_MUSchae_notAac <- c(13, 16, 15, 10, 8, 8, 11, 10, 9, 10, 4, 10, 6, 5, 5, 6, 6, 4, 5, 5, 3, 7, 4, 4)
D_MUSchae_neuDac <- c(11, 6, 8, 9, 9, 9, 9, 9, 7, 8, 8, 6, 4, 5, 7, 7, 3, 5, 8, 3, 4, 4, 6, 3, 5, 3, 1, 1)
D_MUSac_notM_P_A <- c(5, 4, 7, 5, 5, 8, 8, 5, 6, 9, 5, 8, 5, 5, 10, 6, 5, 4, 6, 6, 4, 5, 4, 8, 3, 4, 4, 2, 4, 3, 3, 4, 2, 2, 1, 2)
D_acicula <- c(9, 5, 7, 5, 5, 5, 3, 4, 3, 7, 4, 5)
D_MUSlong_V <- c(14, 14, 13, 11, 10, 5, 9, 7, 11, 9, 8, 7, 11, 11, 10, 8, 4, 13, 10, 9, 5, 4, 9, 7, 7, 10, 9, 8, 6, 9, 9, 8, 7, 5, 4, 5, 8, 6, 6, 5, 9, 6, 6, 6, 6, 4, 4, 6, 6, 6, 6, 5, 5, 6, 6, 5, 4, 5, 4, 4, 4, 3, 5, 4, 4, 4, 3, 3, 3, 3, 4, 3, 3, 2, 3, 2, 2, 2, 2, 1)
D_MUSlong_D <- c(10, 15, 13, 10, 10, 11, 11, 8, 12, 9, 9, 10, 13, 12, 9, 5, 5, 13, 10, 7, 9, 9, 6, 5, 9, 8, 6, 5, 5, 4, 6, 10, 9, 8, 5, 4, 9, 9, 8, 7, 6, 5, 5, 4, 4, 7, 8, 7, 6, 5, 4, 4, 6, 5, 4, 3, 2, 6, 6, 4, 3, 3, 3, 5, 5, 4, 3, 1, 4, 4, 3, 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1)
D_circumchaetal <- c(3, 5, 4, 2, 3, 2, 3, 2, 4, 7, 3, 3, 3, 4, 4, 4, 3, 5, 4, 3, 5, 4, 4, 3, 4, 5, 2, 2, 3, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 2, 2, 1)  

WD_data <- list(WD_circumacicular, WD_MUSchae_notAac, WD_MUSchae_neuDac, WD_MUSac_notM_P_A, WD_acicula, WD_MUSlong_V, WD_MUSlong_D, WD_circumchaetal)

D_data <- list(D_MUSchae_notAac, D_circumacicular,  D_MUSlong_V, D_MUSchae_neuDac, D_MUSlong_D, D_acicula, D_MUSac_notM_P_A, D_circumchaetal)

mean(D_MUSchae_notAac)
mean(D_circumacicular) 
mean(D_MUSlong_V) 
mean(D_MUSchae_neuDac)
mean(D_MUSlong_D)
mean(D_acicula)
mean(D_MUSac_notM_P_A) 
mean(D_circumchaetal)



# Creating basic Jitter
pdf("Weighted_degrees_plot.pdf", width=6, height=6)
par(mar = c(10, 5, 2, 2), family = "sans", font = 2, cex=1) 
boxplot(WD_data, notch=T, add=F, ann=T, names,outpch = NA, 
        col=hcl.colors(8, "Reds2", alpha = 0.8), x_lab=T, las=2,
        varwidth=FALSE, boxlwd=0.5,
        boxwex=0.3,ylab='weighted degree',
        names = c('circumacicular', 'MUSchae_notAac', 'MUSchae_neuDac', 'MUSac_notM/P/A',
                  'aciculae', 'MUSlong_V', 'MUSlong_D',  'circumchaetal')) 
#Setting 'outpch = NA' avoids plotting outliers
#las=2 rotates the x labels
stripchart(WD_data, 
           vertical = TRUE, method = "jitter", jitter=0.1,
           pch = 1, cex=0.4,col = hcl.colors(1,palette = 'Grays',alpha=0.6), bg = "bisque", 
           add = TRUE,
           at=c(0.7,1.7,2.7,3.7,4.7,5.7,6.7,7.7)) 
dev.off()
# Creating basic Jitter
pdf("Degrees_plot.pdf", width=6, height=6)
par(mar = c(10, 5, 2, 2), family = "sans", font = 2, cex=1) 
boxplot(D_data, notch=T, add=F, ann=T, names,outpch = NA, 
        col=hcl.colors(8, "Reds2", alpha = 0.8), x_lab=T, las=2,
        varwidth=FALSE, boxlwd=0.5,
        boxwex=0.3,ylab='weighted degree',
        names = c('MUSchae_notAac','circumacicular',  'MUSlong_V', 'MUSchae_neuDac', 
                  'MUSlong_D', 'aciculae', 'MUSac_notM/P/A', 'circumchaetal')) 
#Setting 'outpch = NA' avoids plotting outliers
#las=2 rotates the x labels
stripchart(D_data, 
           vertical = TRUE, method = "jitter", jitter=0.1,
           pch = 1, cex=0.4,col = hcl.colors(1,palette = 'Grays',alpha=0.6), bg = "bisque", 
           add = TRUE,
           at=c(0.7,1.7,2.7,3.7,4.7,5.7,6.7,7.7)) 
dev.off()
  