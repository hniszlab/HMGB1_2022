################################################################################
## This script contains the code to produce Extended Data Figures 5I-J
## Mensah & Niskanen et al.
## Disruption of nucleolar phase separation in human genetic disease 2022
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

clean_table_nucleoli <- function(x) {         
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
  
  colnames(x)[6:8] <- c("IntensityMean_DAPI_nucleoli","IntensityMean_RFP_nucleoli","IntensityMean_GFPchannel_nucleoli")
  colnames(x)[9:11] <- c("IntensityStd_DAPI_nucleoli","IntensityStd_RFP_nucleoli","IntensityStd_GFPchannel_nucleoli")
  colnames(x)[12:14] <- c("IntensitySum_DAPI_nucleoli","IntensitySum_RFP_nucleoli","IntensitySum_GFPchannel_nucleoli")
  
  return(x)
}

clean_table_out_nucleoli <- function(x) {         
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
  
  colnames(x)[6:8] <- c("IntensityMean_DAPI_out_nucleoli","IntensityMean_RFP_out_nucleoli","IntensityMean_GFPchannel_out_nucleoli")
  colnames(x)[9:11] <- c("IntensityStd_DAPI_out_nucleoli","IntensityStd_RFP_out_nucleoli","IntensityStd_GFPchannel_out_nucleoli")
  colnames(x)[12:14] <- c("IntensitySum_DAPI_out_nucleoli","IntensitySum_RFP_out_nucleoli","IntensitySum_GFPchannel_out_nucleoli")
  
  return(x)
}

nucleis <- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_Nucleis.csv")
nucleis <- clean_table_nuclei(nucleis)

nucleis2 <- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_images2_Nucleis.csv")
nucleis2 <- clean_table_nuclei(nucleis2)

nucleolis<- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_Nucleolis.csv")
nucleolis <- clean_table_nucleoli(nucleolis)

nucleolis2<- read.csv2(file="analysis_files/20210830_cyt_nuc_nuli_run1_images2_Nucleolis.csv")
nucleolis2 <- clean_table_nucleoli(nucleolis2)

out_nucleolis <- read.csv2(file="analysis_files/20210830_cyt_nuc_outNuli_run1_out_nucleolis.csv")
out_nucleolis <- clean_table_out_nucleoli(out_nucleolis)

out_nucleolis2 <- read.csv2(file="analysis_files/20210830_cyt_nuc_outNuli_run1_images2_out_nucleolis.csv")
out_nucleolis2 <- clean_table_out_nucleoli(out_nucleolis2)


nuclei_IF <- rbind(nucleis, nucleis2)
nucleoli_IF <- rbind(nucleolis, nucleolis2)
out_nucleoli_IF <- rbind(out_nucleolis, out_nucleolis2)

comb_table_IF <- cbind(nuclei_IF[,c(1,6:14)], nucleoli_IF[,c(6:14)], out_nucleoli_IF[,c(6:14)])

# Calculate ratio of nucleolar vs non-nucleolar RFP signal (NPM1)
comb_table_IF$nuc_in_out_ratio <- comb_table_IF$IntensityMean_RFP_nucleoli / comb_table_IF$IntensityMean_RFP_out_nucleoli


# Subset by WT or MUT 
comb_table_IF_WT_FL <- comb_table_IF[grep("WT_FL", comb_table_IF$ImageDocumentName..Image.Name), ]
comb_table_IF_WT_FL$protein <- "WT_FL"

comb_table_IF_MUT_FL <- comb_table_IF[grep("MUT_FL", comb_table_IF$ImageDocumentName..Image.Name), ]
comb_table_IF_MUT_FL$protein <- "MUT_FL"

# Combine to one df, no filtering for GFP signal
df_comb_nonfilt <- rbind(comb_table_IF_WT_FL, comb_table_IF_MUT_FL)

# N of nuclei:
nrow(comb_table_IF_WT_FL)
# [1] 114

nrow(comb_table_IF_MUT_FL)
# [1] 123


# Correlation between nucleolar GFP intensity and ratio of RFP (in vs out nucleoli)


cor.test(comb_table_IF_MUT_FL$IntensityMean_GFPchannel_nucleoli, comb_table_IF_MUT_FL$nuc_in_out_ratio)
# Pearson's product-moment correlation
# 
# data:  comb_table_IF_MUT_FL$IntensityMean_GFPchannel_nucleoli and comb_table_IF_MUT_FL$nuc_in_out_ratio
# t = -10.86, df = 119, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7851331 -0.6029415
# sample estimates:
#        cor 
# -0.7055091 

cor.test(comb_table_IF_WT_FL$IntensityMean_GFPchannel_nucleoli, comb_table_IF_WT_FL$nuc_in_out_ratio)
#Pearson's product-moment correlation

#data:  comb_table_IF_WT_FL$IntensityMean_GFPchannel_nucleoli and comb_table_IF_WT_FL$nuc_in_out_ratio
#t = -0.02048, df = 111, p-value = 0.9837
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.1866067  0.1828516
#sample estimates:
#         cor 
# -0.001943897 


# Plot data as scatter plot

ggplot(df_comb_nonfilt, aes(x=IntensityMean_GFPchannel_nucleoli, y=nuc_in_out_ratio, color=protein)) +
  geom_point(shape=16, alpha=0.8) + 
  scale_colour_manual(values = c("#FF0000","#000000")) +
  theme_classic() +
  geom_line(stat="smooth",method = "lm", alpha = 0.6) +
  labs(y ="NPM1 in nucleolus / NPM1 outside nucleolus", x = "Nucleolar GFP intensity") +
  annotate(geom="text", x=120, y=9, label="cor = 0.002, pval = 0.9837", color="black") +
  annotate(geom="text", x=150, y=2.5, label="cor = -0.7055, pval < 2.22e-16", color="red") +
  ylim(0,20)

file <- "Extended_Data_Figure_5I.pdf"
ggsave(file, width = 15 , height = 15, units = "cm", useDingbats =F)



# Correlation between GFP inside nucleoli vs NPM1 outside nucleoli

cor.test(comb_table_IF_MUT_FL$IntensityMean_GFPchannel_nucleoli, comb_table_IF_MUT_FL$IntensityMean_RFP_out_nucleoli)
# Pearson's product-moment correlation
# 
# data:  comb_table_IF_MUT_FL$IntensityMean_GFPchannel_nucleoli and comb_table_IF_MUT_FL$IntensityMean_RFP_out_nucleoli
# t = 6.2393, df = 119, p-value = 6.975e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3489081 0.6200329
# sample estimates:
#       cor 
# 0.4964841 

cor.test(comb_table_IF_WT_FL$IntensityMean_GFPchannel_nucleoli, comb_table_IF_WT_FL$IntensityMean_RFP_out_nucleoli)
# Pearson's product-moment correlation

# data:  comb_table_IF_WT_FL$IntensityMean_GFPchannel_nucleoli and comb_table_IF_WT_FL$IntensityMean_RFP_out_nucleoli
# t = -2.3516, df = 111, p-value = 0.02046
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.38699733 -0.03449944
# sample estimates:
#        cor 
# -0.2178409 

# Plot data as scatter plot

# With lines drawn
ggplot(df_comb_nonfilt, aes(x=IntensityMean_GFPchannel_nucleoli, y=IntensityMean_RFP_out_nucleoli, color=protein)) +
  geom_point(shape=16, alpha=0.8) + 
  scale_colour_manual(values = c("#FF0000","#000000")) +
  theme_classic() +
  geom_line(stat="smooth",method = "lm", alpha = 0.6) +
  labs(y ="NPM1 intensity outside nucleoli", x = "IntensityMean_GFPchannel_nucleoli") +
  annotate(geom="text", x=120, y=0.5, label="cor = -0.2178, pval = 0.02", color="black") +
  annotate(geom="text", x=150, y=6, label="cor = 0.4965, pval = 6.975e-09", color="red") +
  ylim(0,8)

file <- "Extended_Data_Figure_5J.pdf"
ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

