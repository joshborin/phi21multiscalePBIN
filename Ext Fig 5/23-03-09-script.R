library(ggplot2)
library(ggthemes)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggpubr)

setwd('ClinicalTraining/Manuscript/Ext Fig 5')
replay_dat <- read.csv(file = "2022-04-05_phage_titers.csv")
replay_dat <- subset(replay_dat, Type == 'Repeat')
# COLORS ####
colores <- c('#BBBBBB','#FFA525','#68BBFF','#B21F00')
names(colores) <- c('All','OmpC','LamB','OmpF')

# CALC Prop Receptor use ####
# with respect to WT measured on same day
# make index and empty columns
replay_dat$index <- 1:nrow(replay_dat)
replay_dat$PropWT <- 0
# subset for each Day and Population
#  take WT titer for that day
#  calculate PropWT and append to rawdat$PropWT by index
for (i in unique(replay_dat$Timepoint)) {
  for (j in unique(replay_dat$Population)) {
    for (k in unique(replay_dat$Type)) {
      for (l in unique(replay_dat$Rep)) {
        # can set both vals and subset by both at same time
        sub <- subset(replay_dat, Timepoint == i & Population == j & Type == k & Rep == l)
        titer_WT <- sub[sub$Plated_Host == 'WT',]$Titer
        sub$PropWT <- sub$Titer/titer_WT
        replay_dat$PropWT[sub$index] <- sub$PropWT
      }
    }
  }
}
# log10 transform receptor values
replay_dat$logPropWT <- log10(replay_dat$PropWT)
replay_dat$logPropWT[replay_dat$logPropWT == -Inf] <- -6 

# Plots all together ####
recep_plot <- ggplot(data = na.omit(subset(replay_dat, Receptor != 'All')), 
                     aes(x=Timepoint, y=logPropWT, color=Receptor,
                         group=interaction(Receptor, Rep))) +
  geom_line(size = 0.5) + geom_point(size = 1) +
  theme_classic(base_size = 16) +
  labs(y='Frequency of Receptor Use', x='Day') +
  scale_x_continuous(limits=c(0,24), breaks=seq(0,24,3)) +
  scale_y_continuous(limits=c(-6,1), breaks=c(0,-1,-2,-3,-4,-5,-6),
                     labels=c('1','0.1','0.01','0.001','0.0001','0.00001','0')) +
  scale_color_manual(values=colores, breaks=c('LamB','OmpC','OmpF')) +
  facet_wrap(~Rep)
recep_plot

# Plot loop ####
plot_list <- list()
ticker <- 0
for (i in unique(replay_dat$Rep)) {
  rep_sub <- subset(replay_dat, Rep == i)
  ticker <- ticker + 1
  plot_list[[i]] <- ggplot(data = na.omit(subset(rep_sub, Receptor != 'All')), 
                       aes(x=Timepoint, y=logPropWT, color=Receptor)) +
    geom_line(size = 1) +
    theme_classic(base_size = 16) +
    labs(y='Freq. of Receptor Use', x='Day') +
    scale_x_continuous(limits=c(0,24), breaks=seq(0,24,3)) +
    scale_y_continuous(limits=c(-6,1), breaks=c(0,-1,-2,-3,-4,-5,-6),
                       labels=c('1','0.1','0.01','0.001','0.0001','0.00001','0')) +
    scale_color_manual(values=colores, breaks=c('LamB','OmpC','OmpF'))
}

plotmat <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
          plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],
          legend = 'right', common.legend = TRUE,
          ncol = 3, nrow = 3)
ggsave(filename = 'replay-plot.pdf', plotmat,
       width = 12, height = 7)


