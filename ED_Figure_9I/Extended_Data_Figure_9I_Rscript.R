################################################################################
## This script contains the code to produce Extended Data Figure 9I
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(dplyr)

clean_table_nuclei <- function(x) {         
  x <- x[2:nrow(x),]
  
  for (i in c(4:ncol(x))){
    x[,i] <- as.numeric(gsub(",", ".", x[,i]))
  }
  
  colnames(x)[5:7] <- c("IntensityMean_blue_nuclei","IntensityMean_red_nuclei","IntensityMean_green_nuclei")
  colnames(x)[8:10] <- c("IntensityStd_blue_nuclei","IntensityStd_red_nuclei","IntensityStd_green_nuclei")
  colnames(x)[11:13] <- c("IntensitySum_blue_nuclei","IntensitySum_red_nuclei","IntensitySum_green_nuclei")
  
  return(x)
}

nuclei <- read.csv2(file="ED_Figure_9I_data/220919_cyt_nuc_nuli_lowRFP_erode8_Nuclei.csv")
nuclei <- clean_table_nuclei(nuclei)

# Sort out nuclei with low green signal:
nuclei_GFP_filt <- nuclei %>% filter(IntensityMean_green_nuclei > 1)

nuclei_GFP_filt$imagename_nucleus_ID <- paste(nuclei_GFP_filt$ImageDocumentName..Image.Name, nuclei_GFP_filt$ID..ID..I)


# No RFP filter at this point
#nuclei_GFP_filt_RFPabove2 <- nuclei_GFP_filt %>% filter(IntensityMean_red_nuclei >= 2)


clean_table_nucleolis <- function(x) {         
  x <- x[2:nrow(x),]
  
  for (i in c(4:ncol(x))){
    x[,i] <- as.numeric(gsub(",", ".", x[,i]))
  }
  
  colnames(x)[6:8] <- c("IntensityMean_blue_nucleolis","IntensityMean_red_nucleolis","IntensityMean_green_nucleolis")
  colnames(x)[9:11] <- c("IntensityStd_blue_nucleolis","IntensityStd_red_nucleolis","IntensityStd_green_nucleolis")
  colnames(x)[12:14] <- c("IntensitySum_blue_nucleolis","IntensitySum_red_nucleolis","IntensitySum_green_nucleolis")
  
  return(x)
}

clean_table_out_nucleolis <- function(x) {         
  x <- x[2:nrow(x),]
  
  for (i in c(4:ncol(x))){
    x[,i] <- as.numeric(gsub(",", ".", x[,i]))
  }
  
  colnames(x)[6:8] <- c("IntensityMean_blue_out_nucleolis","IntensityMean_red_out_nucleolis","IntensityMean_green_out_nucleolis")
  colnames(x)[9:11] <- c("IntensityStd_blue_out_nucleolis","IntensityStd_red_out_nucleolis","IntensityStd_green_out_nucleolis")
  colnames(x)[12:14] <- c("IntensitySum_blue_out_nucleolis","IntensitySum_red_out_nucleolis","IntensitySum_green_out_nucleolis")
  
  return(x)
}


nucleolis <- read.csv2(file="ED_Figure_9I_data/220919_cyt_nuc_nuli_lowRFP_erode8_Nucleolis.csv")
nucleolis <- clean_table_nucleolis(nucleolis)

out_nucleolis <- read.csv2(file="ED_Figure_9I_data/220919_cyt_nuc_outNuli_lowRFP_erode8_out_nucleolis.csv")
out_nucleolis <- clean_table_out_nucleolis(out_nucleolis)

nucleolis$imagename_nucleus_ID <- paste(nucleolis$ImageDocumentName..Image.Name, nucleolis$ParentID..ID.of.the.parent..I)
out_nucleolis$imagename_nucleus_ID <- paste(out_nucleolis$ImageDocumentName..Image.Name, out_nucleolis$ParentID..ID.of.the.parent..I)


# Select values based on nucleus IDs - include GFP filt
nucleolis1 <- nucleolis %>% filter(imagename_nucleus_ID %in% nuclei_GFP_filt$imagename_nucleus_ID)

out_nucleolis1 <- out_nucleolis %>% filter(imagename_nucleus_ID %in% nuclei_GFP_filt$imagename_nucleus_ID)

# combine tables:
#nucleolis_combined <- rbind(nucleolis1, nucleolis2)
#out_nucleolis_combined <- rbind(out_nucleolis1, out_nucleolis2)

comb_table <- cbind(nucleolis1, out_nucleolis1)


# Remove NaN:s
comb_table_filt <- comb_table[complete.cases(comb_table), ]
comb_table_filt2 <- comb_table_filt
comb_table_filt2$nucleolar_enrichment <- comb_table_filt2$IntensityMean_green_nucleolis/comb_table_filt2$IntensityMean_green_out_nucleolis
comb_table_filt2$log2_nucleolar_enrichment <- log2(comb_table_filt2$nucleolar_enrichment)


df_DVL1_wt <- comb_table_filt2[grep("DVL1_wt", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_DVL1_mut <- comb_table_filt2[grep("DVL1_mut", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_DVL1_WT <- comb_table_filt2[grep("DVL1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_DVL1_MUT <- comb_table_filt2[grep("DVL1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

# Combine 2 names together
df_DVL1_WT <- rbind(df_DVL1_wt, df_DVL1_WT)
df_DVL1_MUT <- rbind(df_DVL1_mut, df_DVL1_MUT)

df_DVL1_WT$sample <- "01_DVL1_wt"
df_DVL1_MUT$sample <- "02_DVL1_mut"

comb_table_filt3 <- rbind(df_DVL1_WT, df_DVL1_MUT)

comb_table_filt3$log2_nucleolar_enrichment <- log2(comb_table_filt3$nucleolar_enrichment)



# Plot values
df_to_plot <- comb_table_filt3[,c(31,32,33)]

p <- ggplot(df_to_plot, aes(x=sample, y=log2_nucleolar_enrichment))

p +  geom_jitter(aes(colour=sample),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  #ylim(0, 30) + 
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  scale_colour_manual(values = c("#A9A9A9","#D22B2B")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(legend.position = "none")
NULL


file <- "Extended_Data_Figure_9I.pdf"
ggsave(file, width = 5 , height = 12, units = "cm")

# N of nuclei (filtered):
nrow(df_DVL1_WT)
#51

nrow(df_DVL1_MUT)
#56

# Print table of results:
print_table <- comb_table_filt3[,c(32,31,33)]

write.table(print_table, file ="DVL1_print_nucleolar_enrichment_quant.txt", row.names = F)


# statistical tests:

t.test(df_DVL1_WT$log2_nucleolar_enrichment, df_DVL1_MUT$log2_nucleolar_enrichment)
# t = -5.4953, df = 81.678, p-value = 4.323e-07



