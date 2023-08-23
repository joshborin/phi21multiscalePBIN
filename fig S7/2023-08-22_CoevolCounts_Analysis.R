library(ggplot2)
library(ggthemes)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggpubr)

setwd('ClinicalTraining/Manuscript/fig S7')
old_dat <- read.csv(file = "prev-P2T0-T21-data.csv")
new_dat <- read.csv(file = "2022-04-05_phage_titers.csv")

# CONTINUE EXPT ####
# Titer plot ####
dat_cont <- subset(new_dat, Type == 'Continue')
dat_cont_p2 <- subset(dat_cont, Timepoint < 25 & Rep == 1)
plot_cont_titer_p2 <- ggplot(subset(dat_cont_p2, Plated_Host == 'WT' & Titer != 'NA'),
                          aes(x = Timepoint,
                              y = log10(as.numeric(Titer)))) +
  geom_line(size = 1) + geom_point(size = 1.4) +
  geom_line(data = old_dat, aes(x = Timepoint,
                                y = log10(as.numeric(Titer_WT))),
            size = 1) +
  geom_point(data = old_dat, aes(x = Timepoint,
                                 y = log10(as.numeric(Titer_WT))),
             size = 1.4) +
  theme_classic(base_size = 16) +
  labs(y = bquote('Phage Titer (log' [10] ~ 'PFU/mL)'),
       x = ('Day')) +
  scale_y_continuous(breaks = 0:11, limits = c(0,11)) +
  scale_x_continuous(breaks = seq(0,24,3), limits = c(0,24)) +
  # hline is the limit of detection
  #geom_hline(yintercept = 2.7, linetype = 'dashed') +
  geom_vline(xintercept = 21, linetype = 'dotted')
plot_cont_titer_p2
ggsave(filename = "cont-titer-plot_pop2.pdf", 
       plot = plot_cont_titer_p2, width = 6, height = 4, bg = 'transparent')





