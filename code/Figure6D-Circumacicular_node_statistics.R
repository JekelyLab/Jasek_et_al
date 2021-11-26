#weighted degrees (WD) and degrees (D) of nodes in the desmosomal connectome
#plot for Figure 6 panel D of Jasek et al 2021 Desmosomal connectome paper
#Gaspar Jekely May 2021

#load packages and functions
library(ggplot2)
library(tidyverse)

#cbind_dif function will return data.frame and fill in NAs to make all columns of the same length
cbind_dif <- function(x = list()){
        # Find max length
        max_length <- max(unlist(lapply(x, length)))
        # Set length of each vector as
        res <- lapply(x, function(x){
                length(x) <- max_length
                return(x)
        })
        return(as.data.frame(res))
}

setwd("/working_directory/")


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


#create tibble with cbind_dif function that will fill in NAs to make all columns of the same length
WD_cbind <- cbind_dif(list(WD_circumacicular, WD_MUSchae_notAac, WD_MUSchae_neuDac, WD_MUSac_notM_P_A, WD_acicula, WD_MUSlong_V, WD_MUSlong_D, WD_circumchaetal))
names(WD_cbind) <- c('circumacicular', 'MUSchae_notAac', 'MUSchae_neuDac', 'MUSac_notM/P/A', 'aciculae', 'MUSlong_V', 'MUSlong_D',  'circumchaetal')
WD_cbind <- tibble(WD_cbind)

#pivot the tibble
WD_cbind <- WD_cbind %>% 
        pivot_longer(everything(),names_to = "cell_type", values_to = "weighted_degree")

#encode cell type vector as a factor
WD_cbind$cell_type <- factor(WD_cbind$cell_type)

#same for unweighted data
#create tibble with cbind_dif function that will fill in NAs to make all columns of the same length
D_cbind <- cbind_dif(list(D_MUSchae_notAac, D_circumacicular,  D_MUSlong_V, D_MUSchae_neuDac, D_MUSlong_D, D_acicula, D_MUSac_notM_P_A, D_circumchaetal))
names(D_cbind) <- c('MUSchae_notAac','circumacicular',  'MUSlong_V', 'MUSchae_neuDac', 'MUSlong_D', 'aciculae', 'MUSac_notM/P/A', 'circumchaetal')
D_cbind <- tibble(D_cbind)
D_cbind

#pivot the tibble
D_cbind <- D_cbind %>% 
        pivot_longer(everything(),names_to = "cell_type", values_to = "weighted_degree")

#encode cell type vector as a factor
D_cbind$cell_type <- factor(D_cbind$cell_type)

#plotting
WD_cbind %>%
        ggplot(mapping=aes(fct_reorder(cell_type, weighted_degree, .fun='mean', na.rm = TRUE, .desc = TRUE), weighted_degree))+
        geom_boxplot(notch=TRUE, fill=hcl.colors(8, palette='Reds2', alpha = 0.6), outlier.shape = NA)+
        geom_jitter(alpha=0.4,width=0.25, height = 0.1, color='grey20', shape=16, size=1.4)+
         #       geom_jitter(mapping=aes(x=(`cell type`),y=`weighted degree`), width=0.2, height = 0.01)+
        theme(panel.background = element_rect(fill = "white", color = "grey"))+
        theme(legend.position="none")+
        ylab("weighted degree") + 
        xlab("")+
        theme(axis.text.x = element_text(angle = 90, family='sans', size=15, face='plain', hjust=1, vjust = 0.5))+
        theme(axis.text.y=element_text(family='sans',size=15))+
        theme(axis.title.y = element_text(size=20))
ggsave('WD_weights_cells.pdf', width=5, height=6)
ggsave('WD_weights_cells.png', width=5, height=6)

#plotting
D_cbind %>%
        ggplot(mapping=aes(fct_reorder(cell_type, weighted_degree, .fun='mean', na.rm = TRUE, .desc = TRUE), weighted_degree))+
        geom_boxplot(notch=TRUE, fill=hcl.colors(8, palette='Reds2', alpha = 0.6), outlier.shape = NA)+
        geom_jitter(alpha=0.4,width=0.25, height = 0.1, color='grey20', shape=16, size=1.4)+
        #       geom_jitter(mapping=aes(x=(`cell type`),y=`weighted degree`), width=0.2, height = 0.01)+
        theme(panel.background = element_rect(fill = "white", color = "grey"))+
        theme(legend.position="none")+
        ylab("weighted degree") + 
        xlab("")+
        theme(axis.text.x = element_text(angle = 90, family='sans', size=15, face='plain', hjust=1, vjust = 0.5))+
        theme(axis.text.y=element_text(family='sans',size=15))+
        theme(axis.title.y = element_text(size=20))
ggsave('D_weights_cells.pdf', width=5, height=6)
ggsave('D_weights_cells.png', width=5, height=6)
