
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(car)

setwd("ClinicalTraining/Manuscript/Figure 2")
dat_propphage <- read.csv(file = "22-02-16-freeze-comm-iso.csv")
dat_propbact <- read.csv(file = "22-12-15-pop-mut-data.csv")

# Manually specify colors for each Treatment
colores <- c('#BBBBBB','#FFA525','#68BBFF','#B21F00')
names(colores) <- c('All','OmpC','LamB','OmpF')

# Phage prop plot ####
dat_propphage$logProp <- log10(dat_propphage$PropWT)
dat_propphage$logProp[dat_propphage$logProp == -Inf] <- -3 

propplot_p <- ggplot(data = subset(dat_propphage, Receptor != 'All'), 
                  aes(x = Timepoint, y = logProp, color = Receptor)) +
  geom_point(size = 2.5) + geom_line(size = 1.5) +
  theme_classic(base_size = 18) +
  labs(y='Freq. of Receptor Use', x='Day') +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,3)) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0),
                     labels = c(0, 0.01, 0.1, 1)) +
  scale_color_manual(values = colores, breaks = c('LamB','OmpC','OmpF'))
propplot_p
ggsave(filename = 'PanelA.pdf', propplot_p, width = 7, height = 4)


# Bacteria prop plot ####
sum_propbact <- dat_propbact %>%
  group_by(Timepoint, Category) %>%
  summarise(
    Freq = sum(Frequency))

propplot_b <- ggplot(sum_propbact, 
               aes(x = Timepoint, y = Freq, color = Category)) +
  geom_point(size = 2.5) + geom_line(size = 1.5) + 
  theme_classic(base_size = 18) +
  labs(y='Freq. of Mutations', x='Day', color = "Receptor") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,3)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_color_manual(values = colores, breaks = c('LamB','OmpC','OmpF'))
propplot_b
ggsave(filename = 'PanelB.pdf', propplot_b,
       width = 7, height = 4)

PanelsAB <- ggarrange(propplot_p, propplot_b, ncol = 1, align = "hv",
                      common.legend = T, legend = "right")
ggsave(filename = 'PanelAB.pdf', PanelsAB, width = 7, height = 7)

# Combined phage/bacteria prop plot ####
ylim_prim <- c(-3.1,0.2)
ylim_sec <- c(-0.03,1.0)

b <- diff(ylim_prim)/diff(ylim_sec)
a <- (ylim_prim[1] - b*ylim_sec[1]) # there was a bug here

propplot_pb <- ggplot(data = subset(dat_propphage, Receptor != 'All'), 
                     aes(x = Timepoint, y = logProp, color = Receptor)) +
  geom_point(size = 2.5) + geom_line(size = 1.5) +
  geom_line(data = sum_propbact,
            aes(x=Timepoint, y=a+(Freq*b),color=Category),
            linetype="dotted",size=1.5) +
  theme_classic(base_size = 18) +
  labs(y='Freq. of Receptor Use', x='Day') +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,3)) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0),
                     labels = c(0, 0.01, 0.1, 1),
                     sec.axis = sec_axis(~ (. - a)/b, name = "Freq. of Receptor Mutations")) +
  scale_color_manual(values = colores, breaks = c('LamB','OmpC','OmpF'))
propplot_pb
ggsave(filename = 'PanelAB-v2.pdf', propplot_pb, width = 8, height = 4)

# This plot is just to get linetype legend
legend_plot <- ggplot(data = subset(dat_propphage, Receptor != 'All'), 
                     aes(x = Timepoint, y = logProp, linetype = Receptor)) +
  geom_point(size = 2.5) + geom_line(size = 1.5) +
  theme_classic(base_size = 18) +
  labs(y='Freq. of Receptor Use', x='Day') +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,3)) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0),
                     labels = c(0, 0.01, 0.1, 1)) +
  theme(legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(2, 'cm'))
legend_plot
ggsave(filename = 'PanelA-legend.pdf', legend_plot, width = 7, height = 4)

# Host range data ####
dat_hostrange <- read.csv(file = "2022-02-24-P2-hostrange.csv")
dat_hostrange <- subset(dat_hostrange, Method == "Spot")
dat_hostrange <- na.omit(dat_hostrange)
dat_hostrange <- subset(dat_hostrange, Timepoint < 15)

colores2 <- c("#808080ff","#F0E442")
names(colores2) <- c("0","1")

# Receptor use data ####
dat_recepuse <- read.csv(file = "2022-03-04-receptor-use.csv")
dat_recepuse <- na.omit(dat_recepuse)

dat_hostrange$ReceptorUse <- dat_recepuse$ReceptorUse[match(dat_hostrange$ID, dat_recepuse$ID)]
dat_hostrange$ID <- str_remove(dat_hostrange$ID, "P2")

# Reorder
dat_hostrange <- dat_hostrange %>%
  mutate(HR.Receptor = fct_relevel(HR.Receptor, 
                                   "OmpF","OmpC","LamB","All"))

level_info <- unique(dat_hostrange$ID)
level_info <- c("Phi21_Anc",
                "T3_WT_1","T3_WT_2","T3_WT_3","T3_C-F-_1","T3_C-F-_2","T3_C-F-_3",
                "T3_L-F-_1","T3_L-F-_2","T3_L-F-_3",
                "T6_WT_1","T6_WT_2","T6_WT_3","T6_C-F-_1","T6_C-F-_2","T6_C-F-_3" ,
                "T6_L-F-_1","T6_L-F-_2","T6_L-F-_3",
                "T9_WT_1","T9_WT_2","T9_WT_3","T9_C-F-_1","T9_C-F-_2","T9_C-F-_3",
                "T9_L-F-_1","T9_L-F-_2","T9_L-F-_3","T9_L-C-_1","T9_L-C-_2","T9_L-C-_3",
                "T12_WT_1","T12_WT_2","T12_WT_3","T12_C-F-_1","T12_C-F-_2","T12_C-F-_3",
                "T12_L-F-_1","T12_L-F-_2","T12_L-F-_3","T12_L-C-_1","T12_L-C-_2","T12_L-C-_3")
dat_hostrange <- dat_hostrange %>%
  mutate(ID = fct_relevel(ID, as.character(level_info)))

# Receptor use plot ####
dat_hostrange$Status[dat_hostrange$Titer > 0] <- 1
dat_hostrange$Status[dat_hostrange$Titer == 0] <- 0

heatmap <- ggplot(subset(dat_hostrange, HR.Receptor != 'All'),
                   aes(y=HR.Receptor, x=ID, fill=as.factor(Status))) + 
  geom_tile(color = 'white', size = 0.2) +
  scale_fill_manual(values = colores2) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Phage Isolate") + ylab("Receptor Use") +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.ticks=element_line(size=0.4),
        axis.text.x.top = element_text(color = 'black', size = 6,
                                       angle=90, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(color = 'black')) +
  geom_hline(yintercept = c(1.5,2.5), color = 'white', size = 1) +
  geom_vline(xintercept = c(1.5,10.5,19.5,31.5,41.5,53.5,62.5,74.5,
                            83.5,91.5,97.5), color = 'white', size = 1)
heatmap
ggsave('host-recep-matrix_v2.pdf', heatmap, height = 2, width = 7)


# Arms Race data ####
dat_ARD <- read.csv(file = "22-10-12-data.csv")
dat_ARD <- subset(dat_ARD, Expt != 'A')

# Ordering
bact_levels <- c('WT','T3_1','JB2','T9_1','JB191')
phage_levels <- c('Anc','P2T3_L-F-_2','P2T9_WT_3')

dat_ARD <- dat_ARD %>%
  mutate(Phage = fct_relevel(Phage,
                             as.character(phage_levels)))
dat_ARD <- dat_ARD %>%
  mutate(PlatedHost = fct_relevel(PlatedHost,
                                  as.character(bact_levels)))

# output binary ####
dat_ARD$Status <- 0
dat_ARD$Status[dat_ARD$Titer > 0] <- 1

# Arms Race plot
ARDplot <- ggplot(dat_ARD,
                  aes(y=PlatedHost,x=Phage,fill=as.factor(Status))) + 
  geom_tile(color = 'white', size = 0.3) +
  scale_fill_manual(values = colores2) +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  xlab("Phage Isolate") + ylab("Bacterial Host") +
  theme(panel.background = element_blank(),
        legend.position = 'right') +
  theme(axis.ticks=element_line(size=0.4),
        axis.text.x.top = element_text(color="black",angle=90, vjust = 0.5, hjust = 0.0),
        axis.text.y = element_text(color="black"),
        axis.title = element_text(color="black",size = 12)) +
  geom_vline(xintercept = c(1.5,2.5), color = 'white', size = 1.5) +
  geom_hline(yintercept = c(2.5,4.5), color = 'white', size = 1.5)
ARDplot
ggsave('PBIN-arms-race.pdf', ARDplot, height = 4, width = 4.5)







