################################################################################
## This script contains the code to produce values Extended Data Figure 5A-B
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(dplyr)

# Function to tidy the table
clean_table_nuclei <- function(x) {
  x <- x[2:nrow(x),]

  x[,6] <- as.numeric(gsub(",", ".", x[,6]))
  x[,7] <- as.numeric(gsub(",", ".", x[,7]))
  x[,8] <- as.numeric(gsub(",", ".", x[,8]))
  x[,9] <- as.numeric(gsub(",", ".", x[,9]))
  x[,10] <- as.numeric(gsub(",", ".", x[,10]))
  x[,11] <- as.numeric(gsub(",", ".", x[,11]))
  x[,12] <- as.numeric(gsub(",", ".", x[,12]))
  x[,13] <- as.numeric(gsub(",", ".", x[,13]))
  x[,14] <- as.numeric(gsub(",", ".", x[,14]))

  colnames(x)[6:8] <- c("IntensityMean_DAPI_nuclei","IntensityMean_RFP_nuclei","IntensityMean_GFPchannel_nuclei")
  colnames(x)[9:11] <- c("IntensityStd_DAPI_nuclei","IntensityStd_RFP_nuclei","IntensityStd_GFPchannel_nuclei")
  colnames(x)[12:14] <- c("IntensitySum_DAPI_nuclei","IntensitySum_RFP_nuclei","IntensitySum_GFPchannel_nuclei")

  return(x)
}

clean_table_cyto <- function(x) {
  x <- x[2:nrow(x),]

  x[,6] <- as.numeric(gsub(",", ".", x[,6]))
  x[,7] <- as.numeric(gsub(",", ".", x[,7]))
  x[,8] <- as.numeric(gsub(",", ".", x[,8]))
  x[,9] <- as.numeric(gsub(",", ".", x[,9]))
  x[,10] <- as.numeric(gsub(",", ".", x[,10]))
  x[,11] <- as.numeric(gsub(",", ".", x[,11]))
  x[,12] <- as.numeric(gsub(",", ".", x[,12]))
  x[,13] <- as.numeric(gsub(",", ".", x[,13]))
  x[,14] <- as.numeric(gsub(",", ".", x[,14]))

  colnames(x)[6:8] <- c("IntensityMean_DAPI_cyto","IntensityMean_RFP_cyto","IntensityMean_GFPchannel_cyto")
  colnames(x)[9:11] <- c("IntensityStd_DAPI_cyto","IntensityStd_RFP_cyto","IntensityStd_GFPchannel_cyto")
  colnames(x)[12:14] <- c("IntensitySum_DAPI_cyto","IntensitySum_RFP_cyto","IntensitySum_GFPchannel_cyto")

  return(x)
}




nuclei1 <- read.csv2(file="analysis_files/set1/Nucleis.csv")
nuclei1 <- clean_table_nuclei(nuclei1)

nuclei2 <- read.csv2(file="analysis_files/set2/Nucleis.csv")
nuclei2 <- clean_table_nuclei(nuclei2)

nuclei3 <- read.csv2(file="analysis_files/set3/Nucleis.csv")
nuclei3 <- clean_table_nuclei(nuclei3)


cyto1 <- read.csv2(file="analysis_files/set1/Cytos.csv")
cyto1 <- clean_table_cyto(cyto1)

cyto2 <- read.csv2(file="analysis_files/set2/Cytos.csv")
cyto2 <- clean_table_cyto(cyto2)

cyto3 <- read.csv2(file="analysis_files/set3/Cytos.csv")
cyto3 <- clean_table_cyto(cyto3)


nuclei_table <- rbind(nuclei1, nuclei2, nuclei3)
cyto_table <- rbind(cyto1, cyto2, cyto3)



comb_table <- cbind(nuclei_table[,c(1,6:14)], cyto_table[,c(6:14)])


comb_table$nuclear_GFP_enrichment <- comb_table$IntensityMean_GFPchannel_nuclei / comb_table$IntensityMean_GFPchannel_cyto

######### normalizing Std by mean intensity
comb_table$normalized_STD <- comb_table$IntensityStd_GFPchannel_nuclei / comb_table$IntensityMean_GFPchannel_nuclei

# Subset by protein
df_WT_FL <- comb_table[grep("WT_FL", comb_table$ImageDocumentName..Image.Name), ]
df_MUT_FL <- comb_table[grep("MUT_FL", comb_table$ImageDocumentName..Image.Name), ]
df_WT_IDR <- comb_table[grep("WT_IDR", comb_table$ImageDocumentName..Image.Name), ]
df_MUT_IDR <- comb_table[grep("MUT_IDR", comb_table$ImageDocumentName..Image.Name), ]

# Add column that states protein name
df_WT_FL$protein <- "01_WT_FL"
df_MUT_FL$protein <- "03_MUT_FL"
df_WT_IDR$protein <- "02_WT_IDR"
df_MUT_IDR$protein <- "04_MUT_IDR"

#Combine back to single df
comb_table <- rbind(df_WT_FL, df_MUT_FL, df_WT_IDR, df_MUT_IDR)


################################################################################
# filter out nuclei with mean EGFP intensity < 5
comb_table_GFP_filt <- comb_table %>% filter(IntensityMean_GFPchannel_nuclei > 5)




################################################################################
# Plot nuclear enrichment

p <- ggplot(comb_table_GFP_filt, aes(x=protein, y=nuclear_GFP_enrichment))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  ylim(0, 30) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c("#A9A9A9","#A9A9A9","#A9A9A9", "#A9A9A9")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL

file <- "Extended_Data_Figure_5A.pdf"
ggsave(file, width = 10 , height = 15, units = "cm")


################################################################################
# Plot normalized standard deviation

p <- ggplot(comb_table_GFP_filt, aes(x=protein, y=normalized_STD))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c("#A9A9A9","#A9A9A9","#A9A9A9", "#A9A9A9")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL


file <- "Extended_Data_Figure_5B.pdf"
ggsave(file, width = 10 , height = 15, units = "cm")


df_WT_FL_GFP_filt <- comb_table_GFP_filt[grep("WT_FL", comb_table_GFP_filt$ImageDocumentName..Image.Name), ]
df_MUT_FL_GFP_filt <- comb_table_GFP_filt[grep("MUT_FL", comb_table_GFP_filt$ImageDocumentName..Image.Name), ]
df_WT_IDR_GFP_filt <- comb_table_GFP_filt[grep("WT_IDR", comb_table_GFP_filt$ImageDocumentName..Image.Name), ]
df_MUT_IDR_GFP_filt <- comb_table_GFP_filt[grep("MUT_IDR", comb_table_GFP_filt$ImageDocumentName..Image.Name), ]


nrow(df_WT_FL_GFP_filt)
# 121
nrow(df_MUT_FL_GFP_filt)
# 102
nrow(df_WT_IDR_GFP_filt)
# 32
nrow(df_MUT_IDR_GFP_filt)
# 38


################################################################################
# T-tests:

# T tests for normalized Std
t.test(df_WT_FL_GFP_filt$normalized_STD, df_WT_IDR$normalized_STD)
# Welch Two Sample t-test
#
# data:  df_WT_FL_GFP_filt$normalized_STD and df_WT_IDR$normalized_STD
# t = -1.7328, df = 46.008, p-value = 0.08983
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.073153403  0.005470387
# sample estimates:
#   mean of x mean of y
# 0.2183397 0.2521812


t.test(df_WT_FL_GFP_filt$normalized_STD, df_MUT_FL$normalized_STD)
# Welch Two Sample t-test
#
# data:  df_WT_FL_GFP_filt$normalized_STD and df_MUT_FL$normalized_STD
# t = -28.056, df = 124.17, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.8492464 -0.7373187
# sample estimates:
#   mean of x mean of y
# 0.2183397 1.0116223

t.test(df_WT_FL_GFP_filt$normalized_STD, df_MUT_IDR$normalized_STD)
# Welch Two Sample t-test
#
# data:  df_WT_FL_GFP_filt$normalized_STD and df_MUT_IDR$normalized_STD
# t = -26.316, df = 45.321, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7067228 -0.6062542
# sample estimates:
#   mean of x mean of y
# 0.2183397 0.8748282

t.test(df_WT_IDR_GFP_filt$normalized_STD, df_MUT_IDR_GFP_filt$normalized_STD)
# Welch Two Sample t-test
#
# data:  df_WT_IDR_GFP_filt$normalized_STD and df_MUT_IDR_GFP_filt$normalized_STD
# t = -38.116, df = 53.278, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7462553 -0.6716508
# sample estimates:
#   mean of x mean of y
# 0.1981026 0.9070556





# T-test for nuclear enrichments:
t.test(df_WT_FL_GFP_filt$nuclear_GFP_enrichment, df_WT_IDR_GFP_filt$nuclear_GFP_enrichment)
# Welch Two Sample t-test
#
# data:  df_WT_FL_GFP_filt$nuclear_GFP_enrichment and df_WT_IDR_GFP_filt$nuclear_GFP_enrichment
# t = 11.526, df = 129.88, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   2.173000 3.073533
# sample estimates:
#   mean of x mean of y
# 3.984104  1.360838


t.test(df_MUT_FL_GFP_filt$nuclear_GFP_enrichment, df_WT_IDR_GFP_filt$nuclear_GFP_enrichment)
# Welch Two Sample t-test
#
# data:  df_MUT_FL_GFP_filt$nuclear_GFP_enrichment and df_WT_IDR_GFP_filt$nuclear_GFP_enrichment
# t = 13.266, df = 104.11, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   4.257230 5.753617
# sample estimates:
#   mean of x mean of y
# 6.366261  1.360838


t.test(df_MUT_IDR_GFP_filt$nuclear_GFP_enrichment, df_WT_IDR_GFP_filt$nuclear_GFP_enrichment)
# Welch Two Sample t-test
#
# data:  df_MUT_IDR_GFP_filt$nuclear_GFP_enrichment and df_WT_IDR_GFP_filt$nuclear_GFP_enrichment
# t = 11.34, df = 37.202, p-value = 1.241e-13
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   8.356243 11.991329
# sample estimates:
#   mean of x mean of y
# 11.534623  1.360838
