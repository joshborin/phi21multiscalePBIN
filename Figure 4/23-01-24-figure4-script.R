# Figure 4

library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ade4)
library(viridis)

setwd("ClinicalTraining/Manuscript/Figure 4")
colores <- c("#440154FF","#FDE725FF")
names(colores) <- c("No","Yes")

# Panel A ####
dattrim <- read.csv(file = "22-08-06-datlong.csv")
dattrim$Phage <- str_replace(dattrim$Phage, "P2","")
dattrim$Status[dattrim$Status == '0'] <- "No"
dattrim$Status[dattrim$Status == '1'] <- "Yes"

# Filter and order ####
bact <- c("T12_12","T15_8","T21_2")
phage <- c("T12_WT_1","T15_WT_1","T18_WT_1")
dattrim <- dattrim[(dattrim$Bacteria %in% bact),]
dattrim <- dattrim[(dattrim$Phage %in% phage),]

dattrim <- dattrim %>%
  mutate(Phage = fct_relevel(Phage,
                             as.character(phage)))
dattrim <- dattrim %>%
  mutate(Bacteria = fct_relevel(Bacteria,
                                as.character(bact)))

# Panel D, Module 1 ####
plot <- ggplot(dattrim, aes(y=Bacteria, x=Phage, fill=Status)) + 
  geom_tile(color = 'white', size = 1) +
  scale_fill_manual(values = colores) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  xlab("Phage Isolates") + ylab("Bacteria Isolates") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.ticks=element_line(size=0.4),
        axis.text.x.top = element_text(size=14,angle=90, vjust = 0.5, hjust = 0.0),
        axis.text.y = element_text(size=14))
plot
ggsave('panelA.pdf', plot, height = 3.5, width = 3)

