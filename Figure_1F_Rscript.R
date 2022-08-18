################################################################################
## This script contains the code to produce Figure 1F
## Mensah & Niskanen et al.
## Disruption of nucleolar phase separation in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(patchwork)

charge_WT <- read.delim(file="EMBOSS_charge_outfile_WT_bin8.txt")
charge_MUT <- read.delim(file="EMBOSS_charge_outfile_MUT_bin8.txt")

p1 <- ggplot(charge_WT, aes(Position, Charge, ymin= 0, ymax = Charge)) + 
        geom_linerange(data = charge_WT, aes(colour = ifelse(Charge <0, "blue", "red")), size=1) +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line()) +
        ylim(-1,1) +
        xlim(0,230) +
        geom_hline(yintercept=0, linetype = "dashed") +
        theme(legend.position="none") + 
        ggtitle("HMGB1 (Wild type)") + 
        theme(plot.title=element_text(hjust=0.5))


p2 <- ggplot(charge_MUT, aes(Position, Charge, ymin= 0, ymax = Charge)) + 
        geom_linerange(data = charge_MUT, aes(colour = ifelse(Charge <0, "blue", "red")), size=1) +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line()) +
        ylim(-1,1) +
        xlim(0,230) +
        geom_hline(yintercept=0, linetype = "dashed") + 
        theme(legend.position="none") +
        ggtitle("HMGB1 (K184Rfs*44)") +
        theme(plot.title=element_text(hjust=0.5))

p1 + p2

file <- "Figure_1F.pdf"
ggsave(file, width = 30 , height = 10, units = "cm")




