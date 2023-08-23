
library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(viridis)
library(ade4)

setwd("ClinicalTraining/Manuscript/fig S3")

# read in each PBIN
orig_dat <- read.csv(file = '22-06-09-matrix-21d.csv')
redo1_dat <- read.csv(file = '22-07-27_data_21d.csv')
redo2_dat <- read.csv(file = '22-09-02_data_21d.csv')
# wide to long
orig_dat_long <- gather(orig_dat, key = Bacteria, value = Status, T3_1:T21_25)
redo1_dat_long <- gather(redo1_dat, key = Bacteria, value = Status, T3_1:T21_23)
redo2_dat_long <- gather(redo2_dat, key = Bacteria, value = Status, T3_1:T21_23)

# MANTEL 1:  Orig vs redo 1 ####
# match for phages and bacteria in common
matchphage <- subset(orig_dat_long, Phage %in% redo1_dat_long$Phage)
matchhost <- subset(matchphage, Bacteria %in% redo1_dat_long$Bacteria)
matchhost[matchhost[,] == 'Y'] <- 1
matchhost[matchhost[,] == 'N'] <- 0
# reorder dataframes to match
new_sorted <- redo1_dat_long[order(redo1_dat_long$Phage, redo1_dat_long$Bacteria),-1]
old_sorted <- matchhost[order(matchhost$Phage, matchhost$Bacteria),]
old_sorted <- old_sorted[,-1]
# spread and make distance matrices
new_wide <- spread(new_sorted,Bacteria,Status)
new_dist <- dist(new_wide)
old_wide <- spread(old_sorted,Bacteria,Status)
old_dist <- dist(old_wide)

mantel.rtest(old_dist, new_dist, nrepet = 999)
# p = 0.001, rho = 0.844


# MANTEL 2: Orig vs redo 2 ####
# match for phages and bacteria in common
matchphage <- subset(orig_dat_long, Phage %in% redo2_dat_long$Phage)
matchhost <- subset(matchphage, Bacteria %in% redo2_dat_long$Bacteria)
matchhost[matchhost[,] == 'Y'] <- 1
matchhost[matchhost[,] == 'N'] <- 0
# reorder dataframes to match
new_sorted <- redo2_dat_long[order(redo2_dat_long$Phage, redo2_dat_long$Bacteria),-1]
old_sorted <- matchhost[order(matchhost$Phage, matchhost$Bacteria),]
old_sorted <- old_sorted[,-1]
# spread and make distance matrices
new_wide <- spread(new_sorted,Bacteria,Status)
new_dist <- dist(new_wide)
old_wide <- spread(old_sorted,Bacteria,Status)
old_dist <- dist(old_wide)

mantel.rtest(old_dist, new_dist, nrepet = 999)
# p = 0.001, rho = 0.832

# PLOT matrices ####
RColorBrewer::brewer.pal(8,'Set1')
colores <- c("#808080ff","#FDE725FF")
names(colores) <- c("No","Yes")
timecolors <- c("#F781BF","#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3","#A65628")

# Panel A ####
# Reorder ####
redo1_dat_long$Bacteria <- str_replace(redo1_dat_long$Bacteria, 'T3_', 'T03_')
redo1_dat_long$Bacteria <- str_replace(redo1_dat_long$Bacteria, 'T6_', 'T06_')
redo1_dat_long$Bacteria <- str_replace(redo1_dat_long$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
redo1_dat_long$Phage <- factor(redo1_dat_long$Phage,
                               levels = unique(redo1_dat_long$Phage))
redo1_dat_long$Bacteria <- factor(redo1_dat_long$Bacteria,
                                  levels = unique(redo1_dat_long$Bacteria))

# HEATMAP ####
redo1_dat_long$Status[redo1_dat_long$Status == '0'] <- "No"
redo1_dat_long$Status[redo1_dat_long$Status == '1'] <- "Yes"

heatmap <- ggplot(redo1_dat_long, aes(y=Bacteria, x=Phage, fill=as.factor(Status))) +
  geom_tile(color = 'white', size = 0.1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  xlab("Phage Isolates") + ylab("Bacterial Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0),
        axis.text.y=element_text(),
        axis.ticks=element_blank()) +
  geom_vline(xintercept = c(9.5,13.5,21.5,31.5,39.5,47.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(7.5,23.5,43.5,50.5,54.5,57.5), color = 'white', size = 1.5)
heatmap
ggsave('PanelA.pdf', heatmap, height = 10, width = 10)


# Panel B ####
# Reorder ####
redo2_dat_long$Bacteria <- str_replace(redo2_dat_long$Bacteria, 'T3_', 'T03_')
redo2_dat_long$Bacteria <- str_replace(redo2_dat_long$Bacteria, 'T6_', 'T06_')
redo2_dat_long$Bacteria <- str_replace(redo2_dat_long$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
redo2_dat_long$Phage <- factor(redo2_dat_long$Phage,
                               levels = unique(redo2_dat_long$Phage))
redo2_dat_long$Bacteria <- factor(redo2_dat_long$Bacteria,
                                  levels = unique(redo2_dat_long$Bacteria))

# HEATMAP ####
redo2_dat_long$Status[redo2_dat_long$Status == '0'] <- "No"
redo2_dat_long$Status[redo2_dat_long$Status == '1'] <- "Yes"

heatmap2 <- ggplot(redo2_dat_long, aes(y=Bacteria, x=Phage, fill=as.factor(Status))) +
  geom_tile(color = 'white', size = 0.1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  xlab("Phage Isolates") + ylab("Bacterial Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0),
        axis.text.y=element_text(),
        axis.ticks=element_blank()) +
  geom_vline(xintercept = c(9.5,13.5,21.5,31.5,39.5,47.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(7.5,22.5,42.5,49.5,53.5,56.5), color = 'white', size = 1.5)
heatmap2
ggsave('PanelB.pdf', heatmap2, height = 10, width = 10)


