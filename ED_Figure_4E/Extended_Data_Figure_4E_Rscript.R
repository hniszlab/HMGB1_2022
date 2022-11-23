################################################################################
## This script contains the code to produce values Extended Data Figure 4E
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(dplyr)

# Droplet analysis from HMGB1 peptide preps
obj <- read.csv2(file = "analysis_files/20220510_B1_peptides_with_RNA_phase_diagram_protein_c_Object.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Remove row # 2 since contains no RNA sample with one object detected
obj <- obj[c(2:nrow(obj)), ]

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

obj$ImageDocumentName..Image.Name

df_PSF <- as.data.frame(cbind(name_vect, vect_a, vect_PSF), stringsAsFactors = F)
colnames(df_PSF) <- c("image_name","area", "PSF")
df_PSF$PSF <- as.numeric(df_PSF$PSF)


# Subset by concentration
df_PSF_01uM <- df_PSF[grep("_0p1uM", df_PSF$image_name), ]
df_PSF_05uM <- df_PSF[grep("_0p5uM", df_PSF$image_name), ]
df_PSF_1uM <- df_PSF[grep("_1uM", df_PSF$image_name), ]
df_PSF_2uM <- df_PSF[grep("_2uM", df_PSF$image_name), ]
df_PSF_5uM <- df_PSF[grep("_5uM", df_PSF$image_name), ]
df_PSF_10uM <- df_PSF[grep("_10uM", df_PSF$image_name), ]

# Add column that states concentration

#df_PSF_01uM$conc <- "01_0.1uM"
df_PSF_05uM$conc <- "02_0.5uM"
df_PSF_1uM$conc <- "03_1uM"
df_PSF_2uM$conc <- "04_2uM"
df_PSF_5uM$conc <- "05_5uM"
df_PSF_10uM$conc <- "06_10uM"

#Combine back to single df
df_PSF <- rbind(df_PSF_05uM, df_PSF_1uM, df_PSF_2uM, df_PSF_5uM, df_PSF_10uM)

# Subset by protein
df_PSF_WT <- df_PSF[grep("WT", df_PSF$image_name), ]
df_PSF_MUT <- df_PSF[grep("MUT", df_PSF$image_name), ]

# Add column that states protein and condition
df_PSF_WT$protein <- "01_WT"
df_PSF_MUT$protein <- "02_MUT"

#Combine back to single df
df_PSF <- rbind(df_PSF_WT, df_PSF_MUT)

# Prepare tables to print
WT_mean_table <- df_PSF_WT %>% group_by(conc) %>% dplyr::summarize(Mean = mean(PSF, na.rm=T), stdev = sd(PSF))
WT_mean_table[is.na(WT_mean_table)] <- 0

MUT_mean_table <- df_PSF_MUT %>% group_by(conc) %>% dplyr::summarize(Mean = mean(PSF, na.rm=T), stdev = sd(PSF))
MUT_mean_table[is.na(MUT_mean_table)] <- 0

WT_mean_table$protein <- c("01_WT")
MUT_mean_table$protein <- c("02_MUT")

# export values
print_table <- rbind(WT_mean_table, MUT_mean_table)

write.table(print_table, file = "Phase_diagram_increasing_protein_conc_mean_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)

# Print all values:
write.table(df_PSF, file = "Phase_diagram_increasing_protein_conc_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)
