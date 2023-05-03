library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(viridis)

setwd("ClinicalTraining/Manuscript/Ext Fig 4")
dat <- read.csv(file = "22-10-05-combined.csv")

# Calc Relative Fitness (W) ####
dat2 <- dat
dat2$index <- seq(1:nrow(dat2))
dat2$W <- NA
for (i in unique(dat2$StrainID)) {
  dat2_s <- subset(dat2, StrainID == i)
  for (j in unique(dat2_s$Rep)) {
    dat2_r <- subset(dat2_s, Rep == j)
    qtiter_t0 <- dat2_r[dat2_r$Time == '0',]$qTiter
    qtiter_t1 <- dat2_r[dat2_r$Time == '1',]$qTiter
    reftiter_t0 <- dat2_r[dat2_r$Time == '0',]$refTiter
    reftiter_t1 <- dat2_r[dat2_r$Time == '1',]$refTiter
    mA <- log(qtiter_t1/qtiter_t0)
    mB <- log(reftiter_t1/reftiter_t0)
    W <- mA/mB
    dat2$W[dat2_r$index] <- W
  }
}

dat2_sub <- subset(dat2, Time == '0')

# Normalize to WT ####
regroup <- dat2_sub
regroup$StrainID[regroup$StrainID == 'T12-11'] <- 'trips'
regroup$StrainID[regroup$StrainID == 'T18-5'] <- 'trips'
regroup$StrainID[regroup$StrainID == 'T21-1'] <- 'trips'
avgdat <- regroup %>%
  group_by(StrainID) %>%
  summarise(meanW = mean(W),
            sdW = sd(W),
            n = n(),
            se = sdW/sqrt(n))
# average relative WT fitness is 1.09 (competitor T6-6)
# so divide all values by 1.092 to normalize to WT


dat2_sub$W <- dat2_sub$W/1.092

# Phylo ordered ####
Wdat <- dat2_sub %>%
  filter(dat2_sub$StrainID != 'WT') %>%
  arrange(W) %>%
  mutate(StrainID = factor(StrainID, levels=c('T3-1','T9-1','T12-11',
                                              'T18-5','T21-1')))

# Statistics ####
anova <- aov(W ~ StrainID, data = subset(Wdat, StrainID != 'WT'))
summary(anova)
TukeyHSD(anova)

W_plot <- ggplot(Wdat, aes(x=StrainID, y=W)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  labs(x= 'Bacterial Isolate',y= 'Relative Fitness (W)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust=0.6, hjust=0.8)) +
  scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1.0)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  annotate(geom = 'text', x = 1, y = 1.07, label = 'A') +
  annotate(geom = 'text', x = 2, y = 1.07, label = 'B') +
  annotate(geom = 'text', x = 3, y = 1.07, label = 'C') +
  annotate(geom = 'text', x = 4, y = 1.07, label = 'C') +
  annotate(geom = 'text', x = 5, y = 1.07, label = 'C')
W_plot
ggsave(W_plot, filename = 'evo-ord-W.pdf', width = 6, height = 4)

# Avg Fitness
regroup <- dat2_sub
regroup$StrainID[regroup$StrainID == 'T12-11'] <- 'trips'
regroup$StrainID[regroup$StrainID == 'T18-5'] <- 'trips'
regroup$StrainID[regroup$StrainID == 'T21-1'] <- 'trips'
avgdat <- regroup %>%
  group_by(StrainID) %>%
  summarise(meanW = mean(W),
            sdW = sd(W),
            n = n(),
            se = sdW/sqrt(n))

