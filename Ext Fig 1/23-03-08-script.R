# Figure 1 Panel A plot script ####

library(ggplot2) 
library(ggpubr)
library(dplyr)
library(tidyverse)
library(ade4)
library(viridis)

setwd("ClinicalTraining/Manuscript/Ext Fig 1")
rawdat <- read.csv(file = '22-06-09-matrix-21d.csv')

# Calculate bacteria resistance range ####
rawdat2 <- rawdat[,-(1:2)]
b_IDs <- colnames(rawdat2)
b_sums <- colSums(rawdat2)
res_dat <- as.data.frame(cbind(b_IDs, b_sums))
res_dat$res_range <- (74 - as.numeric(res_dat$b_sums))

res_dat$b_IDs <- str_replace(res_dat$b_IDs, 'T3_', 'T03_')
res_dat$b_IDs <- str_replace(res_dat$b_IDs, 'T6_', 'T06_')
res_dat$b_IDs <- str_replace(res_dat$b_IDs, 'T9_', 'T09_')

res_dat$Timepoint <- substr(res_dat$b_IDs, 2,3)

b_lm <- lm(res_range~as.numeric(Timepoint), res_dat)
summary(b_lm)
# p < 2e-16

b_anova <- aov(res_range~Timepoint, res_dat)
TukeyHSD(b_anova)

bact_plot <- ggplot(res_dat, aes(x=Timepoint, y=res_range)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0) +
  labs(x= 'Day',y= 'Bacterial Resistance Range') +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits=c(0,80), breaks=c(0,20,40,60,74)) +
  annotate(geom = 'text', x = '03', y = 78, label = 'A') +
  annotate(geom = 'text', x = '06', y = 78, label = 'B') +
  annotate(geom = 'text', x = '09', y = 78, label = 'B') +
  annotate(geom = 'text', x = '12', y = 78, label = 'C') +
  annotate(geom = 'text', x = '15', y = 78, label = 'D') +
  annotate(geom = 'text', x = '18', y = 78, label = 'E') +
  annotate(geom = 'text', x = '21', y = 78, label = 'F')
bact_plot
ggsave(bact_plot, filename = 'PanelA.pdf', width = 4, height = 4)



# Calculate phage host range ####
rawdat$hostsum <- rowSums(rawdat[,-(1:2)])
host_dat <- as.data.frame(cbind(rawdat$Phage, rawdat$hostsum))
colnames(host_dat) <- c('Phage','hostsum')

host_dat$Phage <- str_replace(host_dat$Phage, 'T3_', 'T03_')
host_dat$Phage <- str_replace(host_dat$Phage, 'T6_', 'T06_')
host_dat$Phage <- str_replace(host_dat$Phage, 'T9_', 'T09_')

host_dat$Timepoint <- substr(host_dat$Phage, 4,5)

p_lm <- lm(hostsum~as.numeric(Timepoint), host_dat)
summary(p_lm)
# p < 2e-16

p_anova <- aov(hostsum~Timepoint, host_dat)
summary(p_anova)
TukeyHSD(p_anova)

phage_plot <- ggplot(host_dat, aes(x=Timepoint, y=as.numeric(hostsum))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0) +
  labs(x= 'Day',y= 'Phage Host Range') +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits=c(0,95), breaks=c(0,20,40,60,80,87)) +
  annotate(geom = 'text', x = '03', y = 93, label = 'A') +
  annotate(geom = 'text', x = '06', y = 93, label = 'A') +
  annotate(geom = 'text', x = '09', y = 93, label = 'B') +
  annotate(geom = 'text', x = '12', y = 93, label = 'C') +
  annotate(geom = 'text', x = '15', y = 93, label = 'C') +
  annotate(geom = 'text', x = '18', y = 93, label = 'D') +
  annotate(geom = 'text', x = '21', y = 93, label = 'D')
phage_plot
ggsave(phage_plot, filename = 'PanelB.pdf', width = 4, height = 4)

