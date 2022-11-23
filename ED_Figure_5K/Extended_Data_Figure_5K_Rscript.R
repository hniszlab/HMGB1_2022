################################################################################
## This script contains the code to produce Extended Data Figure 5K
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(dplyr)

# Function to tidy the table
clean_table_nuclei <- function(x) {
  x <- x[2:nrow(x),]

  x[,5] <- as.numeric(gsub(",", ".", x[,5]))
  x[,6] <- as.numeric(gsub(",", ".", x[,6]))
  x[,7] <- as.numeric(gsub(",", ".", x[,7]))

  colnames(x)[5:7] <- c("IntensityMean_GFP_nuclei","IntensityStd_GFP_nuclei","IntensitySum_GFP_nuclei")

  return(x)
}

nuclei1 <- read.csv2(file="data/quant_nuclear_GFP_2_nuc.csv")

nuclei1 <- clean_table_nuclei(nuclei1)

# Normalize Std to mean intensity
nuclei1$STD_norm <- nuclei1$IntensityStd_GFP_nuclei/nuclei1$IntensityMean_GFP_nuclei

# Filter low GFP
nuclei1_GFP_filt <- nuclei1 %>% filter(IntensityMean_GFP_nuclei > 5)

# Subset by protein
HMGB1_WT_FL <- nuclei1_GFP_filt[grep("HMGB1_WT_FL", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_FL <- nuclei1_GFP_filt[grep("HMGB1_MUT_FL", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_WT_delFS <- nuclei1_GFP_filt[grep("HMGB1_WT_delFS", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_Rdel <- nuclei1_GFP_filt[grep("HMGB1_MUT_Rdel", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_RtoA <- nuclei1_GFP_filt[grep("HMGB1_MUT_RtoA", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_RtoK <- nuclei1_GFP_filt[grep("HMGB1_MUT_RtoK", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_KtoR <- nuclei1_GFP_filt[grep("HMGB1_MUT_KtoR", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT_delTail <- nuclei1_GFP_filt[grep("HMGB1_MUT_delTail", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]

# Add column that states protein name
HMGB1_WT_FL$protein <- "01_HMGB1_WT_FL"
HMGB1_MUT_FL$protein <- "02_HMGB1_MUT_FL"
HMGB1_WT_delFS$protein <- "03_HMGB1_WT_delFS"
HMGB1_MUT_Rdel$protein <- "04_HMGB1_MUT_Rdel"
HMGB1_MUT_RtoA$protein <- "05_HMGB1_MUT_RtoA"
HMGB1_MUT_RtoK$protein <- "06_HMGB1_MUT_RtoK"
HMGB1_MUT_KtoR$protein <- "07_HMGB1_MUT_KtoR"
HMGB1_MUT_delTail$protein <- "08_HMGB1_MUT_delTail"


# N of nuclei imaged:
nrow(HMGB1_WT_FL)
# 34
nrow(HMGB1_MUT_FL)
# 17
nrow(HMGB1_WT_delFS)
# 45
nrow(HMGB1_MUT_Rdel)
#23
nrow(HMGB1_MUT_RtoA)
#27
nrow(HMGB1_MUT_RtoK)
#23
nrow(HMGB1_MUT_KtoR)
#13
nrow(HMGB1_MUT_delTail)
#20

#Combine back to single df
comb_table <- rbind(HMGB1_WT_FL, HMGB1_MUT_FL, HMGB1_WT_delFS,
                    HMGB1_MUT_Rdel, HMGB1_MUT_RtoA, HMGB1_MUT_RtoK, HMGB1_MUT_KtoR, HMGB1_MUT_delTail)


# Plot data:

p <- ggplot(comb_table, aes(x=protein, y=STD_norm))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  ylim(0, 3) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c(rep("#A9A9A9",12))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(legend.position = "none", axis.title.x=element_blank())+
  ylab("SD of nuclear EGFP signal normalized by mean intensity") +
  xlab(NA)
  NULL

file <- "Extended_Data_Figure_5K.pdf"
ggsave(file, width = 12 , height = 15, units = "cm")
