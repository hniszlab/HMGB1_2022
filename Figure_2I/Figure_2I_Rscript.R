################################################################################
## This script contains the code to produce Figure 2I
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(patchwork)


################################################################################
# Experiment 1
obj <- read.csv2(file = "analysis_files/peptide_NPM1_exp1_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

unique(obj$ImageDocumentName..Image.Name)
# [1] "20220510_10uM_NPM1_2p5uM_MUT_peptide_image1" "20220510_10uM_NPM1_2p5uM_MUT_peptide_image2" "20220510_10uM_NPM1_2p5uM_MUT_peptide_image3" "20220510_10uM_NPM1_2p5uM_MUT_peptide_image4"
# [5] "20220510_10uM_NPM1_2p5uM_MUT_peptide_image5" "20220510_10uM_NPM1_2p5uM_MUT_peptide_image6" "20220510_10uM_NPM1_2p5uM_MUT_peptide_image7" "20220510_10uM_NPM1_2p5uM_WT_peptide_image1"
# [9] "20220510_10uM_NPM1_2p5uM_WT_peptide_image2"  "20220510_10uM_NPM1_2p5uM_WT_peptide_image3"  "20220510_10uM_NPM1_2p5uM_WT_peptide_image4"  "20220510_10uM_NPM1_2p5uM_WT_peptide_image5"
# [13] "20220510_10uM_NPM1_2p5uM_WT_peptide_image6"  "20220510_10uM_NPM1_5uM_WT_peptide_image1"    "20220510_10uM_NPM1_5uM_WT_peptide_image2"    "20220510_10uM_NPM1_5uM_WT_peptide_image3"
# [17] "20220510_10uM_NPM1_5uM_WT_peptide_image4"    "20220510_10uM_NPM1_5uM_WT_peptide_image5"


# NOTE! "2p5uM" in file name is an error - this was fixed in file names but not in file metadata, that's why an old (incorrect) name is displayed



# Subset by mCh protein:
NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]


R1_df1_wt_0p5 <- NPM1[grep("2p5uM_WT", NPM1$ImageDocumentName..Image.Name), ]
R1_df2_wt_5 <- NPM1[grep("_5uM_WT_", NPM1$ImageDocumentName..Image.Name), ]
R1_df3_mut_0p5 <- NPM1[grep("_2p5uM_MUT", NPM1$ImageDocumentName..Image.Name), ]

R1_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R1_df2_wt_5$protein <- "03_WT_5_uM"
R1_df3_mut_0p5$protein <- "02_MUT_0p5_uM"

comb_table_R1 <- rbind(R1_df1_wt_0p5, R1_df2_wt_5, R1_df3_mut_0p5)


################################################################################
# Experiment #2

obj <- read.csv2(file = "analysis_files/peptide_NPM1_exp2_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))


# Subset by mCh protein:
NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]

R2_df1_wt_0p5 <- NPM1[grep("0p5uM_B1_WT", NPM1$ImageDocumentName..Image.Name), ]
R2_df2_wt_5 <- NPM1[grep("_5uM_B1_WT", NPM1$ImageDocumentName..Image.Name), ]
R2_df3_mut_0p5 <- NPM1[grep("0p5uM_B1_MUT", NPM1$ImageDocumentName..Image.Name), ]
R2_df4_mut_5 <- NPM1[grep("_5uM_B1_MUT", NPM1$ImageDocumentName..Image.Name), ]

R2_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R2_df2_wt_5$protein <- "03_WT_5_uM"
R2_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R2_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R2 <- rbind(R2_df1_wt_0p5, R2_df2_wt_5, R2_df3_mut_0p5, R2_df4_mut_5)


################################################################################
# Experiment #3

obj <- read.csv2(file = "analysis_files/peptide_NPM1_exp3_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]


# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))



# Subset by mCh protein:
NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]

R3_df1_wt_0p5 <- NPM1[grep("WT_0p5uM_", NPM1$ImageDocumentName..Image.Name), ]
R3_df2_wt_5 <- NPM1[grep("WT_5uM_", NPM1$ImageDocumentName..Image.Name), ]
R3_df3_mut_0p5 <- NPM1[grep("MUT_0p5uM_", NPM1$ImageDocumentName..Image.Name), ]
R3_df4_mut_5 <- NPM1[grep("_MUT_5uM_", NPM1$ImageDocumentName..Image.Name), ]

R3_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R3_df2_wt_5$protein <- "03_WT_5_uM"
R3_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R3_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R3 <- rbind(R3_df1_wt_0p5, R3_df2_wt_5, R3_df3_mut_0p5, R3_df4_mut_5)



################################################################################
# Experiment #4           (Considering only NPM1 for now)

obj <- read.csv2(file = "analysis_files/peptide_exp4_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:
NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]

R4_df1_wt_0p5 <- NPM1[grep("0p5uM_B1_WT_", NPM1$ImageDocumentName..Image.Name), ]
R4_df2_wt_5 <- NPM1[grep("_5uM_B1_WT_", NPM1$ImageDocumentName..Image.Name), ]
R4_df3_mut_0p5 <- NPM1[grep("0p5uM_B1_MUT_", NPM1$ImageDocumentName..Image.Name), ]
R4_df4_mut_5 <- NPM1[grep("_5uM_B1_MUT_", NPM1$ImageDocumentName..Image.Name), ]

R4_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R4_df2_wt_5$protein <- "03_WT_5_uM"
R4_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R4_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R4 <- rbind(R4_df1_wt_0p5, R4_df2_wt_5, R4_df3_mut_0p5, R4_df4_mut_5)





################################################################################
# Combine 4 experiments (NPM1)

comb_table_NPM1 <- rbind(comb_table_R1, comb_table_R2, comb_table_R3, comb_table_R4)

p_NPM1 <- ggplot(comb_table_NPM1, aes(x=protein, y=intensity_mean))

p_NPM1 <- p_NPM1 +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 80) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  #geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c(rep(c("#A9A9A9", "#e33232"),2))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("NPM1") +
  xlab("") + ylab("Mean EGFP intensity") +
  theme(legend.position="none")
NULL



t.test(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM", "intensity_mean"],
       comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_NPM1[comb_table_NPM1$protein == "01_WT_0p5_uM", "intensity_mean"] and comb_table_NPM1[comb_table_NPM1$protein == "02_MUT_0p5_uM", "intensity_mean"]
# t = -10.083, df = 740.68, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.9224491 -0.6217812
# sample estimates:
#   mean of x mean of y
# 3.669306  4.441421


t.test(comb_table_NPM1[comb_table_NPM1$protein=="03_WT_5_uM", "intensity_mean"],
       comb_table_NPM1[comb_table_NPM1$protein=="04_MUT_5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_NPM1[comb_table_NPM1$protein == "03_WT_5_uM", "intensity_mean"] and comb_table_NPM1[comb_table_NPM1$protein == "04_MUT_5_uM", "intensity_mean"]
# t = -29.789, df = 414.82, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -11.372509  -9.964524
# sample estimates:
#   mean of x mean of y
# 5.719475 16.387992




# N of observations:

nrow(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM",])
# [1] 350

nrow(comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM",])
# [1] 456

nrow(comb_table_NPM1[comb_table_NPM1$protein=="03_WT_5_uM",])
# [1] 320

nrow(comb_table_NPM1[comb_table_NPM1$protein=="04_MUT_5_uM",])
# [1] 388







################################################################################
# MED1
################################################################################



# MED1 (from experiment 4)

obj <- read.csv2(file = "analysis_files/peptide_exp4_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]


# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:

MED1 <- obj[grep("MED1", obj$ImageDocumentName..Image.Name), ]

R4_df1_wt_0p5 <- MED1[grep("0p5uM_B1_WT_", MED1$ImageDocumentName..Image.Name), ]
R4_df2_wt_5 <- MED1[grep("_5uM_B1_WT_", MED1$ImageDocumentName..Image.Name), ]
R4_df3_mut_0p5 <- MED1[grep("0p5uM_B1_MUT_", MED1$ImageDocumentName..Image.Name), ]
R4_df4_mut_5 <- MED1[grep("_5uM_B1_MUT_", MED1$ImageDocumentName..Image.Name), ]

boxplot(R4_df1_wt_0p5$intensity_mean,
        R4_df2_wt_5$intensity_mean,
        R4_df3_mut_0p5$intensity_mean,
        R4_df4_mut_5$intensity_mean)

R4_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R4_df2_wt_5$protein <- "03_WT_5_uM"
R4_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R4_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R4 <- rbind(R4_df1_wt_0p5, R4_df2_wt_5, R4_df3_mut_0p5, R4_df4_mut_5)




# MED1 (from experiment 5)

obj <- read.csv2(file = "analysis_files/peptide_exp5_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]




# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:

MED1 <- obj[grep("MED1", obj$ImageDocumentName..Image.Name), ]

R5_df1_wt_0p5 <- MED1[grep("0p5uM_B1_WT_", MED1$ImageDocumentName..Image.Name), ]
R5_df2_wt_5 <- MED1[grep("_5uM_B1_WT_", MED1$ImageDocumentName..Image.Name), ]
R5_df3_mut_0p5 <- MED1[grep("0p5uM_B1_MUT_", MED1$ImageDocumentName..Image.Name), ]
R5_df4_mut_5 <- MED1[grep("_5uM_B1_MUT_", MED1$ImageDocumentName..Image.Name), ]

boxplot(R5_df1_wt_0p5$intensity_mean,
        R5_df2_wt_5$intensity_mean,
        R5_df3_mut_0p5$intensity_mean,
        R5_df4_mut_5$intensity_mean)

R5_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R5_df2_wt_5$protein <- "03_WT_5_uM"
R5_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R5_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R5 <- rbind(R5_df1_wt_0p5, R5_df2_wt_5, R5_df3_mut_0p5, R5_df4_mut_5)



# Combine 2 experiments (MED1)

comb_table_MED1 <- rbind(comb_table_R4, comb_table_R5)

p_MED1 <- ggplot(comb_table_MED1, aes(x=protein, y=intensity_mean))

p_MED1 <- p_MED1 +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 150) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  #geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c(rep(c("#A9A9A9", "#e33232"),2))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("MED1") +
  xlab("") + ylab("Mean EGFP intensity") +
  theme(legend.position="none")
NULL




t.test(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM", "intensity_mean"],
       comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM", "intensity_mean"])
# data:  comb_table_MED1[comb_table_MED1$protein == "01_WT_0p5_uM", "intensity_mean"] and comb_table_MED1[comb_table_MED1$protein == "02_MUT_0p5_uM", "intensity_mean"]
# t = 15.902, df = 169.59, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   9.431509 12.105027
# sample estimates:
#   mean of x mean of y
# 14.269815  3.501548


t.test(comb_table_MED1[comb_table_MED1$protein=="03_WT_5_uM", "intensity_mean"],
       comb_table_MED1[comb_table_MED1$protein=="04_MUT_5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_MED1[comb_table_MED1$protein == "03_WT_5_uM", "intensity_mean"] and comb_table_MED1[comb_table_MED1$protein == "04_MUT_5_uM", "intensity_mean"]
# t = 27.009, df = 157.31, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   61.28783 70.95913
# sample estimates:
#   mean of x mean of y
# 78.59928  12.47580







# N of observations:

nrow(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM",])
# [1] 165

nrow(comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM",])
# [1] 133

nrow(comb_table_MED1[comb_table_MED1$protein=="03_WT_5_uM",])
# [1] 153

nrow(comb_table_MED1[comb_table_MED1$protein=="04_MUT_5_uM",])
# [1] 173













################################################################################
# HP1a
################################################################################



# HP1a (from experiment 4)
obj <- read.csv2(file = "analysis_files/peptide_exp4_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:

HP1a <- obj[grep("HP1a", obj$ImageDocumentName..Image.Name), ]

R4_df1_wt_0p5 <- HP1a[grep("0p5uM_B1_WT_", HP1a$ImageDocumentName..Image.Name), ]
R4_df2_wt_5 <- HP1a[grep("_5uM_B1_WT_", HP1a$ImageDocumentName..Image.Name), ]
R4_df3_mut_0p5 <- HP1a[grep("0p5uM_B1_MUT_", HP1a$ImageDocumentName..Image.Name), ]
R4_df4_mut_5 <- HP1a[grep("_5uM_B1_MUT_", HP1a$ImageDocumentName..Image.Name), ]

boxplot(R4_df1_wt_0p5$intensity_mean,
        R4_df2_wt_5$intensity_mean,
        R4_df3_mut_0p5$intensity_mean,
        R4_df4_mut_5$intensity_mean)

R4_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R4_df2_wt_5$protein <- "03_WT_5_uM"
R4_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R4_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R4 <- rbind(R4_df1_wt_0p5, R4_df2_wt_5, R4_df3_mut_0p5, R4_df4_mut_5)



# HP1a (from experiment 5)
obj <- read.csv2(file = "analysis_files/peptide_exp5_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:

HP1a <- obj[grep("HP1a", obj$ImageDocumentName..Image.Name), ]

R5_df1_wt_0p5 <- HP1a[grep("0p5uM_B1_WT_", HP1a$ImageDocumentName..Image.Name), ]
R5_df2_wt_5 <- HP1a[grep("_5uM_B1_WT_", HP1a$ImageDocumentName..Image.Name), ]
R5_df3_mut_0p5 <- HP1a[grep("0p5uM_B1_MUT_", HP1a$ImageDocumentName..Image.Name), ]
R5_df4_mut_5 <- HP1a[grep("_5uM_B1_MUT_", HP1a$ImageDocumentName..Image.Name), ]

boxplot(R5_df1_wt_0p5$intensity_mean,
        R5_df2_wt_5$intensity_mean,
        R5_df3_mut_0p5$intensity_mean,
        R5_df4_mut_5$intensity_mean)

R5_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R5_df2_wt_5$protein <- "03_WT_5_uM"
R5_df3_mut_0p5$protein <- "02_MUT_0p5_uM"
R5_df4_mut_5$protein <- "04_MUT_5_uM"

comb_table_R5 <- rbind(R5_df1_wt_0p5, R5_df2_wt_5, R5_df3_mut_0p5, R5_df4_mut_5)

p <- ggplot(comb_table_R5, aes(x=protein, y=intensity_mean))

p +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  #ylim(0, 3) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  #geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c(rep("#A9A9A9",4))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL






# Combine 2 experiments (HP1a)

comb_table_HP1a <- rbind(comb_table_R4, comb_table_R5)

p_HP1a <- ggplot(comb_table_HP1a, aes(x=protein, y=intensity_mean))

p_HP1a <- p_HP1a +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 80) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  #geom_hline(yintercept=1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c(rep(c("#A9A9A9", "#e33232"),2))) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("HP1a") +
  xlab("") + ylab("Mean EGFP intensity") +
  theme(legend.position="none")
NULL


t.test(comb_table_HP1a[comb_table_HP1a$protein=="01_WT_0p5_uM", "intensity_mean"],
       comb_table_HP1a[comb_table_HP1a$protein=="02_MUT_0p5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_HP1a[comb_table_HP1a$protein == "01_WT_0p5_uM", "intensity_mean"] and comb_table_HP1a[comb_table_HP1a$protein == "02_MUT_0p5_uM", "intensity_mean"]
# t = 6.8395, df = 235.18, p-value = 6.802e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.8368679 1.5140362
# sample estimates:
#   mean of x mean of y
# 5.200994  4.025542


t.test(comb_table_HP1a[comb_table_HP1a$protein=="03_WT_5_uM", "intensity_mean"],
       comb_table_HP1a[comb_table_HP1a$protein=="04_MUT_5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_HP1a[comb_table_HP1a$protein == "03_WT_5_uM", "intensity_mean"] and comb_table_HP1a[comb_table_HP1a$protein == "04_MUT_5_uM", "intensity_mean"]
# t = 3.716, df = 266.79, p-value = 0.0002465
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.045485 3.401868
# sample estimates:
#   mean of x mean of y
# 15.01769  12.79402


# Number of observations:

nrow(comb_table_HP1a[comb_table_HP1a$protein=="01_WT_0p5_uM",])
# [1] 148

nrow(comb_table_HP1a[comb_table_HP1a$protein=="02_MUT_0p5_uM",])
# [1] 176

nrow(comb_table_HP1a[comb_table_HP1a$protein=="03_WT_5_uM",])
# [1] 109

nrow(comb_table_HP1a[comb_table_HP1a$protein=="04_MUT_5_uM",])
# [1] 160




# Plot figure:

p_MED1 + p_NPM1 + p_HP1a



filename <- "Figure_2I.pdf"
ggsave(filename,
       width = 8,
       height = 6,
       dpi = 300,
       useDingbats = FALSE)





################################################################################
# Differences in means:
################################################################################


# NPM1 + 0.5uM peptide
mean(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM","intensity_mean"])
# [1] 3.669306
mean(comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 4.441421

4.441421 / 3.669306
# 1.210425

# NPM1 + 5uM peptide

mean(comb_table_NPM1[comb_table_NPM1$protein=="03_WT_5_uM","intensity_mean"])
# [1] 5.719475
mean(comb_table_NPM1[comb_table_NPM1$protein=="04_MUT_5_uM","intensity_mean"])
# [1] 16.38799

16.38799 / 5.719475
# 2.865296



# MED1 + 0.5uM peptide

mean(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM","intensity_mean"])
#[1] 14.26982
mean(comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 3.501548

3.501548 / 14.26982
# [1]  0.2453814


# MED1 + 5uM peptide

mean(comb_table_MED1[comb_table_MED1$protein=="03_WT_5_uM","intensity_mean"])
# [1] 78.59928
mean(comb_table_MED1[comb_table_MED1$protein=="04_MUT_5_uM","intensity_mean"])
# [1] 12.4758

12.4758 / 78.59928
# 0.1587266



# HP1a + 0.5uM peptide

mean(comb_table_HP1a[comb_table_HP1a$protein=="01_WT_0p5_uM","intensity_mean"])
#[1] 5.200994
mean(comb_table_HP1a[comb_table_HP1a$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 4.025542

4.025542 / 5.200994
# 0.7739947


# HP1a + 5uM peptide

mean(comb_table_HP1a[comb_table_HP1a$protein=="03_WT_5_uM","intensity_mean"])
# [1] 15.01769
mean(comb_table_HP1a[comb_table_HP1a$protein=="04_MUT_5_uM","intensity_mean"])
# [1] 12.79402

12.79402 / 15.01769
# 0.85193
