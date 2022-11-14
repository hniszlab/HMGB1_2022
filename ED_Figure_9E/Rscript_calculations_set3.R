# 2022.09.20

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

nuclei_highRFP <- read.csv2(file="set3_analysis_files/set3_cyt_nuc_nuli_highRFP_Nuclei.csv")
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


nucleolis_highRFP <- read.csv2(file="set3_analysis_files/set3_cyt_nuc_nuli_highRFP_Nucleolis.csv")
nucleolis_lowRFP <- read.csv2(file="set3_analysis_files/set3_cyt_nuc_nuli_lowRFP_Nucleolis.csv")

nucleolis_highRFP <- clean_table_nucleolis(nucleolis_highRFP)
nucleolis_lowRFP <- clean_table_nucleolis(nucleolis_lowRFP)


out_nucleolis_highRFP <- read.csv2(file="set3_analysis_files/set3_cyt_nuc_outNuli_highRFP_out_nucleolis.csv")
out_nucleolis_lowRFP <- read.csv2(file="set3_analysis_files/set3_cyt_nuc_outNuli_lowRFP_out_nucleolis.csv")

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



df_MEN1_WT <- comb_table_filt2[grep("MEN1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_MEN1_MUT <- comb_table_filt2[grep("MEN1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_MEN1_WT$sample <- "01_MEN1_WT"
df_MEN1_MUT$sample <- "02_MEN1_MUT"


df_FOXL2_WT <- comb_table_filt2[grep("FOXL2_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_FOXL2_MUT <- comb_table_filt2[grep("FOXL2_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_FOXL2_WT$sample <- "03_FOXL2_WT"
df_FOXL2_MUT$sample <- "04_FOXL2_MUT"




# NOTE: include SQSTM1 WT DATA FROM EARLIER SET

df_SQSTM1_WT <- comb_table_filt2[grep("20220209_SQSTM1_WT", comb_table_filt2$ImageDocumentName..Image.Name), ]


df_SQSTM1_WT_new <- comb_table_filt2[grep("SQSTM1_WT_new", comb_table_filt2$ImageDocumentName..Image.Name), ]
df_SQSTM1_WT_original <- comb_table_filt2[grep("SQSTM1_WT_original", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_SQSTM1_WT_original2 <- rbind(df_SQSTM1_WT, df_SQSTM1_WT_original) 


df_SQSTM1_MUT <- comb_table_filt2[grep("SQSTM1_MUT", comb_table_filt2$ImageDocumentName..Image.Name), ]

df_SQSTM1_WT_new$sample <- "05_SQSTM1_WT"
df_SQSTM1_WT_original2$sample <- "06_SQSTM1_WT"


df_SQSTM1_MUT$sample <- "07_SQSTM1_MUT"




comb_table_filt3 <- rbind( df_MEN1_WT, df_MEN1_MUT, df_FOXL2_WT, df_FOXL2_MUT, df_SQSTM1_WT_new, df_SQSTM1_WT_original2,df_SQSTM1_MUT)

comb_table_filt3$log2_nucleolar_enrichment <- log2(comb_table_filt3$nucleolar_enrichment)





# Check values by plotting
df_to_plot <- comb_table_filt3[,c(31,32,33)]

ncol(comb_table_filt3)

p <- ggplot(df_to_plot, aes(x=sample, y=log2_nucleolar_enrichment))

p +  geom_jitter(aes(colour=sample),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  #ylim(0, 30) + 
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  #scale_colour_manual(values = rep(c("#A9A9A9","#D22B2B"),12)) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL


# N of nuclei (filtered):


nrow(df_SQSTM1_WT_original2)
# 47
nrow(df_SQSTM1_WT_new)
# 8
nrow(df_SQSTM1_MUT)
# 53

nrow(df_MEN1_WT)
# 63
nrow(df_MEN1_MUT)
# 37

nrow(df_FOXL2_WT)
# 30
nrow(df_FOXL2_MUT)
# 69


# Include only original SQSTM1 data, drop the 8 rows with new version
comb_table_filt4 <- comb_table_filt3[comb_table_filt3$sample != "05_SQSTM1_WT", ]

# Print table of results:
print_table <- comb_table_filt4[,c(32,31,33)]

write.table(print_table, file ="set3_print_nucleolar_enrichment_quant.txt", row.names = F)





# T-tests:
t.test(df_SQSTM1_WT_original2$log2_nucleolar_enrichment, df_SQSTM1_MUT$log2_nucleolar_enrichment)
# t = -3.3963, df = 80.699, p-value = 0.001062

t.test(df_MEN1_WT$log2_nucleolar_enrichment, df_MEN1_MUT$log2_nucleolar_enrichment)
# t = -11.031, df = 90.64, p-value < 2.2e-16

t.test(df_FOXL2_WT$log2_nucleolar_enrichment, df_FOXL2_MUT$log2_nucleolar_enrichment)
# t = -0.25877, df = 74.093, p-value = 0.7965





# Drop SQSTM1 WT from correlation plot analysis:
comb_table_filt5 <- comb_table_filt4[comb_table_filt4$sample != "07_SQSTM1_MUT", ]

# Print table of results:
print_table2 <- comb_table_filt5[,c(32,31,33)]

write.table(print_table2, file ="set3_print_nucleolar_enrichment_quant_noSQSTM1wt.txt", row.names = F)




