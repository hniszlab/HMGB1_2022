################################################################################
## This script contains the code to produce Extended Data Figure 4A
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(patchwork)

WT_mAu <- read.delim(file="WT_mAu_values.txt")

p <- ggplot(data = WT_mAu, aes(x = ml)) +
  geom_line(aes(y = mAU), size = 1) +
  xlim(0,25) +
  xlab("ml") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

WT_plot <- p + geom_segment(aes(x = 12.888525, y = -10, xend = 12.888525, yend = 0)) +
  geom_segment(aes(x = 14.888382, y = -10, xend = 14.888382, yend = 0)) +
  ggtitle("EGFP-HMGB1 WT")

MUT_mAu <- read.delim(file="MUT_mAu_values.txt")

p <- ggplot(data = MUT_mAu, aes(x = ml)) +
  geom_line(aes(y = mAU), size = 1) +
  xlim(0,25) +
  xlab("ml") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

MUT_plot <- p + geom_segment(aes(x = 8.249480, y = -50, xend = 8.249480, yend = 0)) +
  geom_segment(aes(x = 9.749422, y = -50, xend = 9.749422, yend = 0)) +
  ggtitle("EGFP-HMGB1 Mutant")

WT_plot + MUT_plot

file <- "Extended_Data_Figure_4A.pdf"
ggsave(file, width = 25 , height = 8, units = "cm", useDingbats =F)
