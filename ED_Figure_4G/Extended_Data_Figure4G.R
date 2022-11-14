################################################################################
## This script contains the code to produce values Extended Data Figure 4G
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

# Droplet analysis from HMGB1 peptide preps with increasing RNA concentration

library(dplyr)

obj <- read.csv2(file = "analysis_files/20220510_B1_peptides_and_RNA_Object.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset to sample groups:
B1_MUT_samples <- obj[grep("B1_mut", obj$ImageDocumentName..Image.Name), ]

vect_a <- c()
name_vect <- c()
vect_PSF <- c()

# Using total area to calculate fraction of area covered by droplets
# Total area of image:
# 18211810987

# For each image, calculate PSF and print values

for (i in unique(obj$ImageDocumentName..Image.Name)) {

  a <- obj[grep(i, obj$ImageDocumentName..Image.Name), ]
  value <- sum(a$area)
  PSF <- value/(18211810987)

  vect_a <- append(vect_a, value)
  vect_PSF <- append(vect_PSF, PSF)
  name_vect <- append(name_vect, i)

}

df_PSF <- as.data.frame(cbind(name_vect, vect_a, vect_PSF), stringsAsFactors = F)

colnames(df_PSF) <- c("image_name","area", "PSF")

df_PSF$PSF <- as.numeric(df_PSF$PSF)

# Subset by concentration
df_PSF_0ng_RNA <- df_PSF[grep("_0ng_", df_PSF$image_name), ]
df_PSF_2p5ng_RNA <- df_PSF[grep("_2p5ng_", df_PSF$image_name), ]
df_PSF_5ng_RNA <- df_PSF[grep("_5ng_", df_PSF$image_name), ]
df_PSF_10ng_RNA <- df_PSF[grep("_10ng_", df_PSF$image_name), ]
df_PSF_20ng_RNA <- df_PSF[grep("_20ng_", df_PSF$image_name), ]

# Add column that states concentration
df_PSF_0ng_RNA$conc <- "01_0ng"
df_PSF_2p5ng_RNA$conc <- "02_2.5ng"
df_PSF_5ng_RNA$conc <- "03_5ng"
df_PSF_10ng_RNA$conc <- "04_10ng"
df_PSF_20ng_RNA$conc <- "05_20ng"

#Combine back to single df
df_PSF <- rbind(df_PSF_0ng_RNA, df_PSF_2p5ng_RNA, df_PSF_5ng_RNA, df_PSF_10ng_RNA, df_PSF_20ng_RNA)

# Subset by protein
df_PSF_MUT <- df_PSF[grep("mut", df_PSF$image_name), ]

# Add column that states protein and condition
df_PSF_MUT$protein <- "02_MUT"

# Prepare tables to print
MUT_mean_table <- df_PSF_MUT %>% group_by(conc) %>% dplyr::summarize(Mean = mean(PSF, na.rm=T), stdev = sd(PSF))
MUT_mean_table[is.na(MUT_mean_table)] <- 0

MUT_mean_table$protein <- c("02_MUT")

#### Sigmoidal curve fitting:
# export values

write.table(MUT_mean_table, file = "Phase_diagram_increasing_RNA_conc_mean_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)

# Print all values:
write.table(df_PSF, file = "Phase_diagram_increasing_RNA_conc_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)
