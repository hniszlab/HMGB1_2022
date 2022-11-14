################################################################################
## This script contains the code to produce Extended Data Figure 5G
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

nucleis <- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_Nucleis.csv")
nucleis <- clean_table_nuclei(nucleis)

nucleis2 <- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_images2_Nucleis.csv")
nucleis2 <- clean_table_nuclei(nucleis2)


nuclei_IF <- rbind(nucleis, nucleis2)
comb_table_IF <- nuclei_IF[,c(1,6:14)]
comb_table_IF$normalized_SD <- comb_table_IF$IntensityStd_GFPchannel_nuclei / comb_table_IF$IntensityMean_GFPchannel_nuclei

df_WT<- comb_table_IF[grep("WT_FL", comb_table_IF$ImageDocumentName..Image.Name), ]
df_MUT<- comb_table_IF[grep("MUT_FL", comb_table_IF$ImageDocumentName..Image.Name), ]

# Add column that states protein name
df_WT$protein <- "01_WT_FL"
df_MUT$protein <- "02_MUT_FL"


nrow(df_WT)
# 114
nrow(df_MUT)
# 123


# Filter for cells expressing GFP
comb_table_IF_WT_FL_GFP_filt <- df_WT %>% filter(IntensityMean_GFPchannel_nuclei > 5)
comb_table_IF_MUT_FL_GFP_filt <- df_MUT %>% filter(IntensityMean_GFPchannel_nuclei > 5)

nrow(comb_table_IF_WT_FL_GFP_filt)
# 81
nrow(comb_table_IF_MUT_FL_GFP_filt)
# 71

# Test difference in Std
t.test(comb_table_IF_WT_FL_GFP_filt$normalized_SD,
       comb_table_IF_MUT_FL_GFP_filt$normalized_SD)
# data:  comb_table_IF_WT_FL_GFP_filt$normalized_STD and comb_table_IF_MUT_FL_GFP_filt$normalized_STD
# t = -24.839, df = 88.924, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7359616 -0.6269352
# sample estimates:
#   mean of x mean of y
# 0.2276706 0.9091190

#Combine back to single df
comb_table <- rbind(comb_table_IF_WT_FL_GFP_filt, comb_table_IF_MUT_FL_GFP_filt)

p <- ggplot(comb_table, aes(x=protein, y=normalized_SD))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  ylim(0, 2) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c("#A9A9A9","#A9A9A9")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL

file <- "Extended_Data_Figure_5G.pdf"
ggsave(file, width = 10 , height = 15, units = "cm", useDingbats =F)
