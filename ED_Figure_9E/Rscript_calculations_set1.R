# 2021.10.28 

# Quant GFP enrichment at nucleoli vs out_nucleoli


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


nuclei_highRFP <- read.csv2(file="set1_analysis_files/nuc_nuli_highRFP/Nuclei.csv")
nuclei_highRFP <- clean_table_nuclei(nuclei_highRFP)

# Sort out nuclei with low green signal:
nuclei_highRFP_GFP_filt <- nuclei_highRFP %>% filter(IntensityMean_green_nuclei > 5)

nuclei_highRFP_GFP_filt$imagename_nucleus_ID <- paste(nuclei_highRFP_GFP_filt$ImageDocumentName..Image.Name, nuclei_highRFP_GFP_filt$ID..ID..I)

nuclei_highRFP_GFP_filt_RFPabove1 <- nuclei_highRFP_GFP_filt %>% filter(IntensityMean_red_nuclei >= 1)
nuclei_highRFP_GFP_filt_RFPbelow1 <- nuclei_highRFP_GFP_filt %>% filter(IntensityMean_red_nuclei < 1)



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


nucleolis_highRFP <- read.csv2(file="set1_analysis_files/nuc_nuli_highRFP/Nucleolis.csv")
nucleolis_lowRFP <- read.csv2(file="set1_analysis_files/nuc_nuli_lowRFP/Nucleolis.csv")

nucleolis_highRFP <- clean_table_nucleolis(nucleolis_highRFP)
nucleolis_lowRFP <- clean_table_nucleolis(nucleolis_lowRFP)


out_nucleolis_highRFP <- read.csv2(file="set1_analysis_files/nuc_out_nuli_highRFP/out_nucleolis.csv")
out_nucleolis_lowRFP <- read.csv2(file="set1_analysis_files/nuc_out_nuli_lowRFP/out_nucleolis.csv")

out_nucleolis_highRFP <- clean_table_out_nucleolis(out_nucleolis_highRFP)
out_nucleolis_lowRFP <- clean_table_out_nucleolis(out_nucleolis_lowRFP)


nucleolis_highRFP$imagename_nucleus_ID <- paste(nucleolis_highRFP$ImageDocumentName..Image.Name, nucleolis_highRFP$ParentID..ID.of.the.parent..I)
nucleolis_lowRFP$imagename_nucleus_ID <- paste(nucleolis_lowRFP$ImageDocumentName..Image.Name, nucleolis_lowRFP$ParentID..ID.of.the.parent..I)
out_nucleolis_highRFP$imagename_nucleus_ID <- paste(out_nucleolis_highRFP$ImageDocumentName..Image.Name, out_nucleolis_highRFP$ParentID..ID.of.the.parent..I)
out_nucleolis_lowRFP$imagename_nucleus_ID <- paste(out_nucleolis_lowRFP$ImageDocumentName..Image.Name, out_nucleolis_lowRFP$ParentID..ID.of.the.parent..I)



# Select values based on nucleus IDs

nucleolis1 <- nucleolis_highRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPabove1$imagename_nucleus_ID)
nucleolis2 <- nucleolis_lowRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPbelow1$imagename_nucleus_ID)


out_nucleolis1 <- out_nucleolis_highRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPabove1$imagename_nucleus_ID)
out_nucleolis2 <- out_nucleolis_lowRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPbelow1$imagename_nucleus_ID)


# combine tables:
nucleolis_combined <- rbind(nucleolis1, nucleolis2)
out_nucleolis_combined <- rbind(out_nucleolis1, out_nucleolis2)

comb_table <- cbind(nucleolis_combined, out_nucleolis_combined)


# Remove NaN:s
comb_table_filt <- comb_table[complete.cases(comb_table), ]

# remove overexposed images (mentioned in file name):
comb_table_filt2 <- comb_table_filt[!grepl("overexp", comb_table_filt$ImageDocumentName..Image.Name),]


comb_table_filt2$nucleolar_enrichment <- comb_table_filt2$IntensityMean_green_nucleolis/comb_table_filt2$IntensityMean_green_out_nucleolis
comb_table_filt2$log2_nucleolar_enrichment <- log2(comb_table_filt2$nucleolar_enrichment)



df_HMGB3_WT <- comb_table_filt2[grep("HMGB3_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_HMGB3_MUT <- comb_table_filt2[grep("HMGB3_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_HMGB3_WT$sample <- "01_HMGB3_WT"
df_HMGB3_MUT$sample <- "02_HMGB3_MUT"

df_FOXC1_WT <- comb_table_filt2[grep("FOXC1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_FOXC1_MUT <- comb_table_filt2[grep("FOXC1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_FOXC1_WT$sample <- "03_FOXC1_WT"
df_FOXC1_MUT$sample <- "04_FOXC1_MUT"

df_FOXF1_WT <- comb_table_filt2[grep("FOXF1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_FOXF1_MUT <- comb_table_filt2[grep("FOXF1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_FOXF1_WT$sample <- "05_FOXF1_WT"
df_FOXF1_MUT$sample <- "06_FOXF1_MUT"

df_MYOD_WT <- comb_table_filt2[grep("MYOD_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_MYOD_MUT <- comb_table_filt2[grep("MYOD_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_MYOD_WT$sample <- "07_MYOD_WT"
df_MYOD_MUT$sample <- "08_MYOD_MUT"


df_RAX_WT <- comb_table_filt2[grep("RAX_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_RAX_MUT <- comb_table_filt2[grep("RAX_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]


df_RAX_WT$sample <- "09_RAX_WT"
df_RAX_MUT$sample <- "10_RAX_MUT"

df_RUNX1_WT <- comb_table_filt2[grep("RUNX1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_RUNX1_MUT <- comb_table_filt2[grep("RUNX1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_RUNX1_WT$sample <- "11_RUNX1_WT"
df_RUNX1_MUT$sample <- "12_RUNX1_MUT"



comb_table_filt3 <- rbind(df_HMGB3_WT, df_HMGB3_MUT, df_FOXC1_WT, df_FOXC1_MUT, df_FOXF1_WT, df_FOXF1_MUT, 
                          df_MYOD_WT, df_MYOD_MUT, df_RAX_WT, df_RAX_MUT, df_RUNX1_WT, df_RUNX1_MUT)


comb_table_filt3$log2_nucleolar_enrichment <- log2(comb_table_filt3$nucleolar_enrichment)



# N of nuclei (filtered):
nrow(df_HMGB3_WT)
#24
nrow(df_HMGB3_MUT)
#20

nrow(df_FOXC1_WT)
# 21
nrow(df_FOXC1_MUT)
# 28

nrow(df_FOXF1_WT)
# 17
nrow(df_FOXF1_MUT)
# 17

nrow(df_MYOD_WT)
# 33
nrow(df_MYOD_MUT)
# 30

nrow(df_RAX_WT)
# 31
nrow(df_RAX_MUT)
# 35

nrow(df_RUNX1_WT)
# 28
nrow(df_RUNX1_MUT)
# 31


# Print table of results:
print_table <- comb_table_filt3[,c(32,31,33)]
write.table(print_table, file ="set1_print_nucleolar_enrichment_quant.txt", row.names = F)



# T-tests:

t.test(df_HMGB3_WT$log2_nucleolar_enrichment, df_HMGB3_MUT$log2_nucleolar_enrichment)
#t = -12.32, df = 29.085, p-value = 4.562e-13

t.test(df_FOXC1_WT$log2_nucleolar_enrichment, df_FOXC1_MUT$log2_nucleolar_enrichment)
# t = -20.424, df = 37.388, p-value < 2.2e-16

t.test(df_FOXF1_WT$log2_nucleolar_enrichment, df_FOXF1_MUT$log2_nucleolar_enrichment)
#t = -8.7708, df = 27.412, p-value = 1.922e-09

t.test(df_MYOD_WT$log2_nucleolar_enrichment, df_MYOD_MUT$log2_nucleolar_enrichment)
#t = -7.1594, df = 50.536, p-value = 3.196e-09

t.test(df_RAX_WT$log2_nucleolar_enrichment, df_RAX_MUT$log2_nucleolar_enrichment)
# t = -16.617, df = 50.495, p-value < 2.2e-16

t.test(df_RUNX1_WT$log2_nucleolar_enrichment, df_RUNX1_MUT$log2_nucleolar_enrichment)
#t = -4.036, df = 39.498, p-value = 0.0002422




