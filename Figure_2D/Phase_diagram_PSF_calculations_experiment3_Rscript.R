################################################################################
## This script contains the code to produce values for Figure 2D
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

# Droplet analysis from HMGB1 protein preps with (SEC cleaned preps - experiment 3)
library(dplyr)

obj <- read.csv2(file = "analysis_results_experiment3/20220506_SEC_clean_phase_diagram_3_Object.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset to sample groups:
B1_WT_samples <- obj[grep("B1_WT", obj$ImageDocumentName..Image.Name), ]
B1_MUT_samples <- obj[grep("B1_MUT", obj$ImageDocumentName..Image.Name), ]

df1 <- B1_WT_samples[grep("0p1uM", B1_WT_samples$ImageDocumentName..Image.Name), ]
df2 <- B1_WT_samples[grep("0p5uM", B1_WT_samples$ImageDocumentName..Image.Name), ]
df3 <- B1_WT_samples[grep("_1uM", B1_WT_samples$ImageDocumentName..Image.Name), ]
df4 <- B1_WT_samples[grep("_2uM", B1_WT_samples$ImageDocumentName..Image.Name), ]
df5 <- B1_WT_samples[grep("_5uM", B1_WT_samples$ImageDocumentName..Image.Name), ]
df6 <- B1_WT_samples[grep("_10uM", B1_WT_samples$ImageDocumentName..Image.Name), ]

df1_mut <- B1_MUT_samples[grep("0p1uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]
df2_mut <- B1_MUT_samples[grep("0p5uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]
df3_mut <- B1_MUT_samples[grep("_1uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]
df4_mut <- B1_MUT_samples[grep("_2uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]
df5_mut <- B1_MUT_samples[grep("_5uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]
df6_mut <- B1_MUT_samples[grep("_10uM", B1_MUT_samples$ImageDocumentName..Image.Name), ]

vect_a <- c()
name_vect <- c()
vect_PSF <- c()

# using total area to calculate fraction of droplets within this area

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
colnames(df_PSF) <- c("image_name","sum_of_intensity_sums", "PSF")
df_PSF$PSF <- as.numeric(df_PSF$PSF)

# Subset by concentration
df_PSF_01uM <- df_PSF[grep("_0p1uM", df_PSF$image_name), ]
df_PSF_05uM <- df_PSF[grep("_0p5uM", df_PSF$image_name), ]
df_PSF_1uM <- df_PSF[grep("_1uM", df_PSF$image_name), ]
df_PSF_2uM <- df_PSF[grep("_2uM", df_PSF$image_name), ]
df_PSF_5uM <- df_PSF[grep("_5uM", df_PSF$image_name), ]
df_PSF_10uM <- df_PSF[grep("_10uM", df_PSF$image_name), ]

# Add column that states concentration
df_PSF_01uM$conc <- "01_0.1uM"
df_PSF_05uM$conc <- "02_0.5uM"
df_PSF_1uM$conc <- "03_1uM"
df_PSF_2uM$conc <- "04_2uM"
df_PSF_5uM$conc <- "05_5uM"
df_PSF_10uM$conc <- "06_10uM"

#Combine back to single df
df_PSF <- rbind(df_PSF_01uM, df_PSF_05uM, df_PSF_1uM, df_PSF_2uM, df_PSF_5uM, df_PSF_10uM)

# Subset by protein
df_PSF_WT <- df_PSF[grep("WT", df_PSF$image_name), ]
df_PSF_MUT <- df_PSF[grep("MUT", df_PSF$image_name), ]

# Add column that states protein and condition
df_PSF_WT$protein <- "01_WT"
df_PSF_MUT$protein <- "02_MUT"

#Combine back to single df
df_PSF <- rbind(df_PSF_WT, df_PSF_MUT)




################################################################################
# Prepare tables to print
WT_mean_table <- df_PSF_WT %>% group_by(conc) %>% dplyr::summarize(Mean = mean(PSF, na.rm=T), stdev = sd(PSF))
WT_mean_table[is.na(WT_mean_table)] <- 0

MUT_mean_table <- df_PSF_MUT %>% group_by(conc) %>% dplyr::summarize(Mean = mean(PSF, na.rm=T), stdev = sd(PSF))
MUT_mean_table[is.na(MUT_mean_table)] <- 0

WT_mean_table$protein <- c("01_WT")
MUT_mean_table$protein <- c("02_MUT")


# export values
print_table <- rbind(WT_mean_table, MUT_mean_table)
write.table(print_table, file = "Experiment_3_mean_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)

# Print all values:
write.table(df_PSF, file = "Experiment_3_PSF_values.txt", sep = "\t",
            row.names = F, col.names = TRUE, quote = F)
