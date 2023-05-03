# Figure 1 Panel A plot script ####

library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ade4)
library(viridis)
#library(NatParksPalettes)

setwd("ClinicalTraining/Manuscript/Ext Fig 2")
RColorBrewer::brewer.pal(8,'Set1')
colores <- c("#808080ff","#FDE725FF")
names(colores) <- c("No","Yes")
timecolors <- c("#F781BF","#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3","#A65628")


b_list <- read.csv(file = "strainlist-bacteria.csv")
p_list <- subset(read.csv(file = "strainlist-phage.csv"), Phage != 'P2T12_C-F-_3')

b_colors <- timecolors
names(b_colors) <- levels(b_list$Color)
p_colors <- timecolors
names(p_colors) <- levels(p_list$Color)

# Panel A ####
rawdat <- read.csv(file = "22-06-09_data_21d_long.csv")
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
  geom_vline(xintercept = c(9.5,18.5,30.5,41.5,53.5,62.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(25.5,51.5,76.5,89.5,102.5,115.5), color = 'white', size = 1.5)
heatmap
ggsave('PanelA.pdf', heatmap, height = 16, width = 10)

# Panel B Brim Matrix ####
brimdat <- read.csv("23-03-07_brimdat_long.csv")

brimdat$Bacteria <- str_replace(brimdat$Bacteria, 'T3_', 'T03_')
brimdat$Bacteria <- str_replace(brimdat$Bacteria, 'T6_', 'T06_')
brimdat$Bacteria <- str_replace(brimdat$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
brimdat$Phage <- factor(brimdat$Phage, levels = unique(brimdat$Phage))
brimdat$Bacteria <- factor(brimdat$Bacteria, levels = unique(brimdat$Bacteria))

brim_blist <- data_frame(unique(brimdat$Bacteria))
brim_blist$Color <- 0
colnames(brim_blist)[1] <- 'Bacteria'
brim_blist$Color[match(b_list$Bacteria, brim_blist$Bacteria)] <- b_list$Color

brim_plist <- data_frame(unique(brimdat$Phage))
brim_plist$Color <- 0
colnames(brim_plist)[1] <- 'Phage'
brim_plist$Color[match(p_list$Phage, brim_plist$Phage)] <- p_list$Color

# HEATMAP ####
brimdat$Status[brimdat$Status == '0'] <- "No"
brimdat$Status[brimdat$Status == '1'] <- "Yes"
brim_heatmap <- ggplot(brimdat, aes(y=Bacteria, x=Phage, fill=as.factor(Status))) +
  geom_tile(color = 'white', size = 0.1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  xlab("Phage Isolates") + ylab("Bacterial Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(color = p_colors[brim_plist$Color], angle = 90, hjust = 0),
        axis.text.y=element_text(color = b_colors[rev(brim_blist$Color)]),
        axis.ticks=element_blank())
brim_heatmap
ggsave('PanelB.pdf', brim_heatmap, height = 16, width = 10)

setdiff(brim_plist$Phage, rawdat2$Phage)

