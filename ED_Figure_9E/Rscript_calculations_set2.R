# 2022.02.23 

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

nuclei_highRFP <- read.csv2(file="set2_analysis_files/cyt_nuc_nuli_highRFP_Nuclei.csv")
nuclei_highRFP <- clean_table_nuclei(nuclei_highRFP)

# Sort out nuclei with low green signal:
nuclei_highRFP_GFP_filt <- nuclei_highRFP %>% filter(IntensityMean_green_nuclei > 5)

nuclei_highRFP_GFP_filt$imagename_nucleus_ID <- paste(nuclei_highRFP_GFP_filt$ImageDocumentName..Image.Name, nuclei_highRFP_GFP_filt$ID..ID..I)

nuclei_highRFP_GFP_filt_RFPabove8 <- nuclei_highRFP_GFP_filt %>% filter(IntensityMean_red_nuclei >= 8)
nuclei_highRFP_GFP_filt_RFPbelow8 <- nuclei_highRFP_GFP_filt %>% filter(IntensityMean_red_nuclei < 8)



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


nucleolis_highRFP <- read.csv2(file="set2_analysis_files/cyt_nuc_nuli_highRFP_Nucleolis.csv")
nucleolis_lowRFP <- read.csv2(file="set2_analysis_files/cyt_nuc_nuli_lowRFP_Nucleolis.csv")

nucleolis_highRFP <- clean_table_nucleolis(nucleolis_highRFP)
nucleolis_lowRFP <- clean_table_nucleolis(nucleolis_lowRFP)


out_nucleolis_highRFP <- read.csv2(file="set2_analysis_files/cyt_nuc_Out_nuli_highRFP_out_nucleolis.csv")
out_nucleolis_lowRFP <- read.csv2(file="set2_analysis_files/cyt_nuc_Out_nuli_lowRFP_out_nucleolis.csv")

out_nucleolis_highRFP <- clean_table_out_nucleolis(out_nucleolis_highRFP)
out_nucleolis_lowRFP <- clean_table_out_nucleolis(out_nucleolis_lowRFP)


nucleolis_highRFP$imagename_nucleus_ID <- paste(nucleolis_highRFP$ImageDocumentName..Image.Name, nucleolis_highRFP$ParentID..ID.of.the.parent..I)
nucleolis_lowRFP$imagename_nucleus_ID <- paste(nucleolis_lowRFP$ImageDocumentName..Image.Name, nucleolis_lowRFP$ParentID..ID.of.the.parent..I)
out_nucleolis_highRFP$imagename_nucleus_ID <- paste(out_nucleolis_highRFP$ImageDocumentName..Image.Name, out_nucleolis_highRFP$ParentID..ID.of.the.parent..I)
out_nucleolis_lowRFP$imagename_nucleus_ID <- paste(out_nucleolis_lowRFP$ImageDocumentName..Image.Name, out_nucleolis_lowRFP$ParentID..ID.of.the.parent..I)



# Select values based on nucleus IDs

nucleolis1 <- nucleolis_highRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPabove8$imagename_nucleus_ID)
nucleolis2 <- nucleolis_lowRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPbelow8$imagename_nucleus_ID)


out_nucleolis1 <- out_nucleolis_highRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPabove8$imagename_nucleus_ID)
out_nucleolis2 <- out_nucleolis_lowRFP %>% filter(imagename_nucleus_ID %in% nuclei_highRFP_GFP_filt_RFPbelow8$imagename_nucleus_ID)


# combine tables:
nucleolis_combined <- rbind(nucleolis1, nucleolis2)
out_nucleolis_combined <- rbind(out_nucleolis1, out_nucleolis2)

comb_table <- cbind(nucleolis_combined, out_nucleolis_combined)


# Remove NaN:s
comb_table_filt <- comb_table[complete.cases(comb_table), ]


comb_table_filt2 <- comb_table_filt

comb_table_filt2$nucleolar_enrichment <- comb_table_filt2$IntensityMean_green_nucleolis/comb_table_filt2$IntensityMean_green_out_nucleolis

comb_table_filt2$log2_nucleolar_enrichment <- log2(comb_table_filt2$nucleolar_enrichment)



df_PHOX2B_WT <- comb_table_filt2[grep("PHOX2B_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_PHOX2B_MUT <- comb_table_filt2[grep("PHOX2B_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_PHOX2B_WT$sample <- "01_PHOX2B_WT"
df_PHOX2B_MUT$sample <- "02_PHOX2B_MUT"

df_SOX2_WT <- comb_table_filt2[grep("SOX2_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_SOX2_MUT <- comb_table_filt2[grep("SOX2_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_SOX2_WT$sample <- "03_SOX2_WT"
df_SOX2_MUT$sample <- "04_SOX2_MUT"


df_CALR_WT <- comb_table_filt2[grep("CALR_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_CALR_MUT <- comb_table_filt2[grep("CALR_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_CALR_WT$sample <- "05_CALR_WT"
df_CALR_MUT$sample <- "06_CALR_MUT"


df_SQSTM1_WT <- comb_table_filt2[grep("SQSTM1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_SQSTM1_MUT <- comb_table_filt2[grep("SQSTM1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_SQSTM1_WT$sample <- "07_SQSTM1_WT"
df_SQSTM1_MUT$sample <- "08_SQSTM1_MUT"


df_MEN1_WT <- comb_table_filt2[grep("MEN1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_MEN1_MUT <- comb_table_filt2[grep("MEN1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_MEN1_WT$sample <- "09_MEN1_WT"
df_MEN1_MUT$sample <- "10_MEN1_MUT"


df_FOXL2_WT <- comb_table_filt2[grep("FOXL2_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_FOXL2_MUT <- comb_table_filt2[grep("FOXL2_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_FOXL2_WT$sample <- "11_FOXL2_WT"
df_FOXL2_MUT$sample <- "12_FOXL2_MUT"


comb_table_filt3 <- rbind(df_PHOX2B_WT, df_PHOX2B_MUT, df_SOX2_WT, df_SOX2_MUT, df_CALR_WT, df_CALR_MUT, 
                          df_SQSTM1_WT, df_SQSTM1_MUT, df_MEN1_WT, df_MEN1_MUT, df_FOXL2_WT, df_FOXL2_MUT)

comb_table_filt3$log2_nucleolar_enrichment <- log2(comb_table_filt3$nucleolar_enrichment)

# N of nuclei (filtered):

nrow(df_PHOX2B_WT)
#41
nrow(df_PHOX2B_MUT)
#33

nrow(df_SOX2_WT)
# 72
nrow(df_SOX2_MUT)
# 57

nrow(df_CALR_WT)
# 83
nrow(df_CALR_MUT)
# 36

nrow(df_SQSTM1_WT)
# 38
nrow(df_SQSTM1_MUT)
# 42

nrow(df_MEN1_WT)
# 46
nrow(df_MEN1_MUT)
# 43

nrow(df_FOXL2_WT)
# 43
nrow(df_FOXL2_MUT)
# 53


# Print table of results:
print_table <- comb_table_filt3[,c(32,31,33)]

write.table(print_table, file ="set2_print_nucleolar_enrichment_quant.txt", row.names = F)




# T-tests:

t.test(df_PHOX2B_WT$log2_nucleolar_enrichment, df_PHOX2B_MUT$log2_nucleolar_enrichment)
# t = -18.878, df = 67.065, p-value < 2.2e-16

t.test(df_SOX2_WT$log2_nucleolar_enrichment, df_SOX2_MUT$log2_nucleolar_enrichment)
# t = -6.9694, df = 92.808, p-value = 4.548e-10

t.test(df_CALR_WT$log2_nucleolar_enrichment, df_CALR_MUT$log2_nucleolar_enrichment)
# t = -5.5613, df = 86.233, p-value = 2.952e-07

t.test(df_SQSTM1_WT$log2_nucleolar_enrichment, df_SQSTM1_MUT$log2_nucleolar_enrichment)
# t = -2.7964, df = 68.644, p-value = 0.006698

t.test(df_MEN1_WT$log2_nucleolar_enrichment, df_MEN1_MUT$log2_nucleolar_enrichment)
# t = -4.9725, df = 50.748, p-value = 7.944e-06

t.test(df_FOXL2_WT$log2_nucleolar_enrichment, df_FOXL2_MUT$log2_nucleolar_enrichment)
# t = 5.0638, df = 93.629, p-value = 2.052e-06




