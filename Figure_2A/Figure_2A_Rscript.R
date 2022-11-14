################################################################################
## This script contains the code to produce Figure 2A
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)

WT_PONDR_result <- read.delim(file="WT_FL_VLXT.txt")
colnames(WT_PONDR_result) <- c("Num","Residue","VLXT","protein")

WT_PONDR_result$protein <- "WT"

MUT_PONDR_result <- read.delim(file="MUT_FL_VLXT.txt")
colnames(MUT_PONDR_result) <- c("Num","Residue","VLXT")
MUT_PONDR_result$protein <- "Mutant"

df <- rbind(WT_PONDR_result, MUT_PONDR_result)

ggplot(data=df, aes(x=Num, y=VLXT, group = protein)) +
  geom_line(aes(color = protein, linetype=protein), size=1.5) +
  labs(y ="PONDR Score (VL-XT)", x = "Amino acid position") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1.2) +
  theme(legend.position="none") +
  ggtitle("HMGB1") +
  # Add HMG box and IDR coordinates
  geom_segment(aes(x = 9, y = 1.15, xend = 79, yend = 1.15), colour = "blue", size = 3) +
  geom_segment(aes(x = 95, y = 1.15, xend = 163, yend = 1.15), colour = "blue", size = 3) +
  geom_segment(aes(x = 135, y = 1.2, xend = 226, yend = 1.2), colour = "orange", size = 3) +
  # Add FS position
  geom_segment(aes(x = 182, y = 1.1, xend = 186, yend = 1.1), colour = "red", size = 3)

file <- "Figure_2A.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)
