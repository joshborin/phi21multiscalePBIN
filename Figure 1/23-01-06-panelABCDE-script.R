# Figure 1 Panel A plot script ####

library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ade4)
library(viridis)

setwd("ClinicalTraining/Manuscript/Figure 1")
rawdat <- read.csv(file = "22-06-09-matrix-21d.csv")
colores <- c("#808080ff","#FDE725FF")
names(colores) <- c("No","Yes")

# Panel A ####
# data formatting ####
# wide to long
rawdat_long <- gather(rawdat, key = Bacteria, value = Status, T3_1:T21_25)
write.csv(rawdat_long, '22-06-09_data_21d_long.csv')

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
  geom_vline(xintercept = c(9.5,18.5,30.5,41.5,53.5,62.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(25.5,51.5,76.5,89.5,102.5,115.5), color = 'white', size = 1.5)
heatmap
ggsave('PanelA-matrix.pdf', heatmap, height = 8, width = 6)

# Panel B Brim Matrix ####
brimdat <- read.csv("22-06-09_d21_LPbrim.csv", check.names = F)
colnames(brimdat)[1] <- "Bacteria"
brim_long <- gather(brimdat, Phage, Status, 2:75, factor_key=TRUE)
write.csv(brim_long, '23-03-07_brimdat_long.csv')

# New Matrix Interaction Network ####
# Reorder ####
brim_long$Bacteria <- str_replace(brim_long$Bacteria, 'T3_', 'T03_')
brim_long$Bacteria <- str_replace(brim_long$Bacteria, 'T6_', 'T06_')
brim_long$Bacteria <- str_replace(brim_long$Bacteria, 'T9_', 'T09_')
# Order phages by order in df
brim_long$Phage <- factor(brim_long$Phage, levels = unique(brim_long$Phage))
brim_long$Bacteria <- factor(brim_long$Bacteria, levels = unique(brim_long$Bacteria))

# HEATMAP ####
brim_long$Status[brim_long$Status == '0'] <- "No"
brim_long$Status[brim_long$Status == '1'] <- "Yes"
brim_heatmap <- ggplot(brim_long,
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
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())
brim_heatmap
ggsave('PanelB-matrix2.pdf', brim_heatmap, height = 8, width = 6)


# Panel C ####
dat <- read.csv("23-04-13-mod-nest-data.csv")

barcols <- c("black","white")
names(colores) <- c("1","2&3")

# Combined nestedness/modularity plot ####
ylim_prim <- c(0,0.8)
ylim_sec <- c(0,0.32)

b <- diff(ylim_prim)/diff(ylim_sec)
a <- (ylim_prim[1] - b*ylim_sec[1]) # there was a bug here

i_plot <- ggplot(data = subset(dat, Metric == 'Nestedness'), 
                    aes(x = Metric, y = Value, fill = Module)) +
  geom_col(size = 1, position = position_dodge(), color = 'black') +
  theme_classic(base_size = 20) +
  theme(legend.position = "right") +
  labs(y='NODF', x='Metric') +
  scale_fill_manual(values = barcols) +
  #scale_x_discrete(breaks = c('A','B'),
  #                 labels = c('Nestedness','Modularity')) +
  scale_y_continuous(limits = c(0,0.8))
i_plot

ii_plot <- ggplot(data = subset(dat, Metric == 'Modularity'), 
                 aes(x = Metric, y = Value, fill = Module)) +
  geom_col(size = 1, position = position_dodge(), color = 'black') +
  theme_classic(base_size = 20) +
  theme(legend.position = "right") +
  labs(y='Qb', x='Metric') +
  scale_fill_manual(values = barcols) +
  #scale_x_discrete(breaks = c('A','B'),
  #                 labels = c('Nestedness','Modularity')) +
  scale_y_continuous(limits = c(0,0.32))
ii_plot

panelC <- ggarrange(i_plot, ii_plot, common.legend = T)
ggsave('panelC.pdf', panelC, height = 5, width = 6)



# Panel D/E ####
dattrim <- read.csv(file = "22-08-06-datlong.csv")
dattrim$Phage <- str_replace(dattrim$Phage, "P2","")
dattrim$Status[dattrim$Status == '0'] <- "No"
dattrim$Status[dattrim$Status == '1'] <- "Yes"

# Filter and order ####
bact <- c("T12_12","T15_9","T15_18","T15_8","T15_1","T15_17","T21_2",
          "T15_14","T15_7","T18_1","T18_22","T18_25","T18_14","T18_11","T18_16","T18_12")
phage <- c("T12_WT_1","T15_WT_1","T15_C-F-_2","T18_WT_1","T18_L-C-_1",
           "T21_L-F-_1","T21_L-F-_2","T21_WT_3","T21_C-F-_2","T21_L-C-_2")
dattrim <- dattrim[(dattrim$Bacteria %in% bact),]
dattrim <- dattrim[(dattrim$Phage %in% phage),]

dattrim <- dattrim %>%
  mutate(Phage = fct_relevel(Phage,
                             as.character(phage)))
dattrim <- dattrim %>%
  mutate(Bacteria = fct_relevel(Bacteria,
                                as.character(bact)))

# Panel D, Module 1 ####
bact_mod1 <- c("T12_12","T15_9","T15_18","T15_8","T15_1","T15_17","T21_2")
phage_mod1 <- c("T12_WT_1","T15_WT_1","T15_C-F-_2","T18_WT_1","T18_L-C-_1")
dat_mod1 <- dattrim[(dattrim$Bacteria %in% bact_mod1),]
dat_mod1 <- dat_mod1[(dat_mod1$Phage %in% phage_mod1),]

plot_mod1 <- ggplot(dat_mod1,
                     aes(y=Bacteria, x=Phage, fill=Status)) + 
  geom_tile(color = 'white', size = 0.1) +
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
plot_mod1
ggsave('module1.pdf', plot_mod1, height = 4, width = 3)


# Panel E, Module 2 ####
bact_mod2 <- c("T15_14","T15_7","T18_1","T18_14","T18_11","T18_16","T18_12")
phage_mod2 <- c("T21_L-F-_1","T21_L-F-_2","T21_WT_3","T21_C-F-_2","T21_L-C-_2")
dat_mod2 <- dattrim[(dattrim$Bacteria %in% bact_mod2),]
dat_mod2 <- dat_mod2[(dat_mod2$Phage %in% phage_mod2),]

plot_mod2 <- ggplot(dat_mod2,
                    aes(y=Bacteria, x=Phage, fill=Status)) + 
  geom_tile(color = 'white', size = 0.1) +
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
plot_mod2
ggsave('module2.pdf', plot_mod2, height = 4, width = 3)


