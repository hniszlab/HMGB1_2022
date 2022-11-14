################################################################################
## This script contains the code to produce Extended Data Figure 6F
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(dplyr)

clean_table_nuclei <- function(x) {
  x <- x[2:nrow(x),]

  x[,5] <- as.numeric(gsub(",", ".", x[,5]))
  x[,6] <- as.numeric(gsub(",", ".", x[,6]))
  x[,7] <- as.numeric(gsub(",", ".", x[,7]))
  x[,8] <- as.numeric(gsub(",", ".", x[,8]))
  x[,9] <- as.numeric(gsub(",", ".", x[,9]))
  x[,10] <- as.numeric(gsub(",", ".", x[,10]))

  colnames(x)[6] <- c("IntensityMean_GFP_nuclei")
  colnames(x)[8] <- c("IntensityStd_GFP_nuclei")
  colnames(x)[10] <- c("IntensitySum_GFP_nuclei")

  return(x)
}

nuclei1 <- read.csv2(file="Nuclei.csv")
nuclei1 <- clean_table_nuclei(nuclei1)

# Normalize Std to mean intensity
nuclei1$STD_norm <- nuclei1$IntensityStd_GFP_nuclei/nuclei1$IntensityMean_GFP_nuclei

# Filter low GFP  (also removes NO DOX samples)
nuclei1_GFP_filt <- nuclei1 %>% filter(IntensityMean_GFP_nuclei > 3)

# Subset by protein
HMGB1_WT <- nuclei1_GFP_filt[grep("WT_DOX", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_MUT <- nuclei1_GFP_filt[grep("MUT_DOX", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]
HMGB1_delTail <- nuclei1_GFP_filt[grep("MUTdeltail", nuclei1_GFP_filt$ImageDocumentName..Image.Name), ]

# Add column that states protein name
HMGB1_WT$protein <- "01_HMGB1_WT_FL"
HMGB1_MUT$protein <- "02_HMGB1_MUT_FL"
HMGB1_delTail$protein <- "03_HMGB1_MUT_Patchless"

# N of nuclei considered:
nrow(HMGB1_WT)
# 94
nrow(HMGB1_MUT)
# 57
nrow(HMGB1_delTail)
# 51

#Conbine back to single df
comb_table <- rbind(HMGB1_WT, HMGB1_MUT, HMGB1_delTail)

# Test difference in normalized SD
t.test(HMGB1_WT$STD_norm, HMGB1_MUT$STD_norm)
#Welch Two Sample t-test

#data:  HMGB1_WT$STD_norm and HMGB1_MUT$STD_norm
#t = -35.11, df = 67.058, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.223837 -1.092174
#sample estimates:
#  mean of x mean of y
#0.3812011 1.5392070

t.test(HMGB1_WT$STD_norm, HMGB1_delTail$STD_norm)
# Welch Two Sample t-test
#
# data:  HMGB1_WT$STD_norm and HMGB1_delTail$STD_norm
# t = -35.198, df = 61.202, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.161305 -1.036456
# sample estimates:
#   mean of x mean of y
# 0.3812011 1.4800813

t.test(HMGB1_MUT$STD_norm, HMGB1_delTail$STD_norm)

# Welch Two Sample t-test
#
# data:  HMGB1_MUT$STD_norm and HMGB1_delTail$STD_norm
# t = 1.3675, df = 106, p-value = 0.1744
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02659535  0.14484664
# sample estimates:
#   mean of x mean of y
# 1.539207  1.480081

# Plot data:
p <- ggplot(comb_table, aes(x=protein, y=STD_norm))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  ylim(0, 3) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  #geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c(rep("#A9A9A9",7))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL

file <- "Extended_Data_Figure_6F.pdf"
ggsave(file, width = 10 , height = 15, units = "cm", useDingbats=FALSE)
