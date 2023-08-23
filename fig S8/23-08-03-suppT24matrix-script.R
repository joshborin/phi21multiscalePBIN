library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ade4)
library(viridis)

setwd("ClinicalTraining/Manuscript/fig S8")
rawdat <- read.csv(file = "22-06-09-matrix-d24.csv")
colores <- c("#808080ff","#FDE725FF")
names(colores) <- c("No","Yes")

colnames(rawdat)[1] <- 'Phage'

# data formatting ####
# wide to long
rawdat_long <- gather(rawdat, key = Bacteria, value = Status, T3_1:T24_24)
write.csv(rawdat_long, '22-06-09_data_24d_long.csv')

# New Matrix Interaction Network ####
# Reorder ####
rawdat_long$Bacteria <- str_replace(rawdat_long$Bacteria, 'T3_', 'T03_')
rawdat_long$Bacteria <- str_replace(rawdat_long$Bacteria, 'T6_', 'T06_')
rawdat_long$Bacteria <- str_replace(rawdat_long$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
rawdat_long$Phage <- factor(rawdat_long$Phage, levels = unique(rawdat_long$Phage))
rawdat_long$Bacteria <- factor(rawdat_long$Bacteria, levels = unique(rawdat_long$Bacteria))

# HEATMAP ####
rawdat_long$Status[rawdat_long$Status == '0'] <- "No"
rawdat_long$Status[rawdat_long$Status == '1'] <- "Yes"
heatmap <- ggplot(rawdat_long,
                  aes(y=Bacteria, x=Phage, fill=as.factor(Status))) +
  geom_tile(color = 'white', size = 0.1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  xlab("Phage Isolates") + ylab("Bacterial Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        #axis.text.x=element_text(angle = 90),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank()) +
  geom_vline(xintercept = c(9.5,18.5,30.5,41.5,53.5,62.5,74.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(24.5,49.5,75.5,100.5,113.5,126.5,139.5), color = 'white', size = 1.5)
heatmap
ggsave('T24-matrix.pdf', heatmap, height = 9, width = 7)



############################
timecolors <- c("#F781BF","#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3","#A65628", "#000000")

b_list <- read.csv(file = "strainlist-bacteria.csv")
p_list <- subset(read.csv(file = "strainlist-phage.csv"), Phage != 'P2T12_C-F-_3')

b_colors <- timecolors
names(b_colors) <- levels(b_list$Color)
p_colors <- timecolors
names(p_colors) <- levels(p_list$Color)

# Panel A ####
rawdat <- read.csv(file = "22-06-09_data_24d_long.csv")
rawdat2 <- rawdat
# Reorder ####
rawdat2$Bacteria <- str_replace(rawdat2$Bacteria, 'T3_', 'T03_')
rawdat2$Bacteria <- str_replace(rawdat2$Bacteria, 'T6_', 'T06_')
rawdat2$Bacteria <- str_replace(rawdat2$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
rawdat2$Phage <- factor(rawdat2$Phage, levels = unique(rawdat2$Phage))
rawdat2$Bacteria <- factor(rawdat2$Bacteria, levels = unique(rawdat2$Bacteria))

# HEATMAP ####
rawdat2$Status[rawdat2$Status == '0'] <- "No"
rawdat2$Status[rawdat2$Status == '1'] <- "Yes"
heatmap <- ggplot(rawdat2, aes(y=Bacteria, x=Phage, fill=as.factor(Status))) +
  geom_tile(color = 'white', size = 0.1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  xlab("Phage Isolates") + ylab("Bacterial Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(color = p_colors[p_list$Color], angle = 90, hjust = 0),
        axis.text.y=element_text(color = b_colors[rev(b_list$Color)]),
        axis.ticks=element_blank()) +
  geom_vline(xintercept = c(9.5,18.5,30.5,41.5,53.5,62.5,74.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(24.5,49.5,75.5,100.5,113.5,126.5,139.5), color = 'white', size = 1.5)
heatmap
ggsave('T24-matrix-labeled.pdf', heatmap, height = 18, width = 12)











