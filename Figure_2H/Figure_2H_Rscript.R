################################################################################
## This script contains the code to produce Figure 2H
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)
library(patchwork)

# Droplet analysis - co-mixing of HMGB1 and mCherry droplets

################################################################################
# Experiment #1

obj <- read.csv2(file = "analysis_files/proteins_experiment_1_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]

# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))

# Subset by mCh protein:  NPM1
NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]

R1_NPM1_df1_wt_0p5 <- NPM1[grep("WT_0p5uM_", NPM1$ImageDocumentName..Image.Name), ]
R1_NPM1_df2_mut_0p5 <- NPM1[grep("MUT_0p5uM_", NPM1$ImageDocumentName..Image.Name), ]

R1_NPM1_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R1_NPM1_df2_mut_0p5$protein <- "02_MUT_0p5_uM"


# Subset by mCh protein: MED1
MED1 <- obj[grep("MED1", obj$ImageDocumentName..Image.Name), ]

R1_MED1_df1_wt_0p5 <-MED1[grep("WT_0p5uM_", MED1$ImageDocumentName..Image.Name), ]
R1_MED1_df2_mut_0p5 <- MED1[grep("MUT_0p5uM_", MED1$ImageDocumentName..Image.Name), ]

R1_MED1_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R1_MED1_df2_mut_0p5$protein <- "02_MUT_0p5_uM"



################################################################################
# Experiment #2

obj <- read.csv2(file = "analysis_files/proteins_experiment_2_mCh_droplet.csv", stringsAsFactors = F)


# rm extra row
obj <- obj[2:nrow(obj),]



# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))


# Subset by mCh protein:

################################################################################
# NPM1

NPM1 <- obj[grep("NPM1", obj$ImageDocumentName..Image.Name), ]

R2_NPM1_df1_wt_0p5 <- NPM1[grep("0p5uM_HMGB1_WT", NPM1$ImageDocumentName..Image.Name), ]
R2_NPM1_df2_mut_0p5 <- NPM1[grep("0p5uM_HMGB1_MUT", NPM1$ImageDocumentName..Image.Name), ]


R2_NPM1_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R2_NPM1_df2_mut_0p5$protein <- "02_MUT_0p5_uM"


################################################################################
# MED1
MED1 <- obj[grep("MED1", obj$ImageDocumentName..Image.Name), ]

R2_MED1_df1_wt_0p5 <- MED1[grep("0p5uM_HMGB1_WT", MED1$ImageDocumentName..Image.Name), ]
R2_MED1_df2_mut_0p5 <- MED1[grep("0p5uM_HMGB1_MUT", MED1$ImageDocumentName..Image.Name), ]

R2_MED1_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R2_MED1_df2_mut_0p5$protein <- "02_MUT_0p5_uM"


################################################################################
# HP1a
HP1a <- obj[grep("HP1a", obj$ImageDocumentName..Image.Name), ]

R2_HP1a_df1_wt_0p5 <- HP1a[grep("0p5uM_HMGB1_WT", HP1a$ImageDocumentName..Image.Name), ]
R2_HP1a_df2_mut_0p5 <- HP1a[grep("0p5uM_HMGB1_MUT", HP1a$ImageDocumentName..Image.Name), ]


R2_HP1a_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R2_HP1a_df2_mut_0p5$protein <- "02_MUT_0p5_uM"





################################################################################
# Experiment #3
obj <- read.csv2(file = "analysis_files/proteins_experiment_3_mCh_droplet.csv", stringsAsFactors = F)

# rm extra row
obj <- obj[2:nrow(obj),]


# Replace , with ., create a new column
obj$intensity_mean <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R))
obj$intensity_sum <- as.numeric(gsub(",", ".", obj$IntensitySum1_Ch1..Intensity.Sum.of.channel..Ch1...R))
obj$intensity_Std <- as.numeric(gsub(",", ".", obj$IntensityStd_Ch1..Intensity.Standard.Deviation.of.channel..Ch1...R))
obj$intensity_mean_mCh <- as.numeric(gsub(",", ".", obj$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R))
obj$area <- as.numeric(gsub(",", ".", obj$Area..Area..R))


# Subset by mCh protein:

################################################################################
# HP1a
HP1a <- obj[grep("HP1a", obj$ImageDocumentName..Image.Name), ]

R3_HP1a_df1_wt_0p5 <- HP1a[grep("0p5uM_B1_WT", HP1a$ImageDocumentName..Image.Name), ]
R3_HP1a_df2_mut_0p5 <- HP1a[grep("0p5uM_B1_MUT", HP1a$ImageDocumentName..Image.Name), ]

R3_HP1a_df1_wt_0p5$protein <- "01_WT_0p5_uM"
R3_HP1a_df2_mut_0p5$protein <- "02_MUT_0p5_uM"

################################################################################







################################################################################
# Combining datasets:
################################################################################


################################################################################
# NPM1:
comb_table_NPM1 <- rbind(R1_NPM1_df1_wt_0p5, R1_NPM1_df2_mut_0p5,
                         R2_NPM1_df1_wt_0p5, R2_NPM1_df2_mut_0p5)


p1 <- ggplot(comb_table_NPM1, aes(x=protein, y=intensity_mean))

p1 <- p1 +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 100) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c("#A9A9A9", "#e33232")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("NPM1") +
  xlab("") + ylab("Mean EGFP intensity") +
  theme(legend.position="none")
NULL


t.test(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM", "intensity_mean"],
       comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM", "intensity_mean"])
# Welch Two Sample t-test
#data:  comb_table_NPM1[comb_table_NPM1$protein == "01_WT_0p5_uM", "intensity_mean"] and comb_table_NPM1[comb_table_NPM1$protein == "02_MUT_0p5_uM", "intensity_mean"]
# t = -8.2296, df = 209.05, p-value = 2.002e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -15.456104  -9.482218
# sample estimates:
#   mean of x mean of y
# 6.701919 19.171080





################################################################################
# MED1

comb_table_MED1 <- rbind(R1_MED1_df1_wt_0p5, R1_MED1_df2_mut_0p5,
                         R2_MED1_df1_wt_0p5, R2_MED1_df2_mut_0p5)


p2 <- ggplot(comb_table_MED1, aes(x=protein, y=intensity_mean))

p2 <- p2 +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 100) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c("#A9A9A9", "#e33232")) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle("MED1") +
  xlab("") + ylab("Mean EGFP intensity") +
  theme(legend.position="none")
NULL



t.test(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM", "intensity_mean"],
       comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM", "intensity_mean"])
# Welch Two Sample t-test
#
# data:  comb_table_MED1[comb_table_MED1$protein == "01_WT_0p5_uM", "intensity_mean"] and comb_table_MED1[comb_table_MED1$protein == "02_MUT_0p5_uM", "intensity_mean"]
# t = -2.9694, df = 220.66, p-value = 0.003313
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -8.768633 -1.772542
# sample estimates:
#   mean of x mean of y
# 21.62651  26.89710



################################################################################
# HP1a

comb_table_HP1a <- rbind(R2_HP1a_df1_wt_0p5, R2_HP1a_df2_mut_0p5,
                         R3_HP1a_df1_wt_0p5, R3_HP1a_df2_mut_0p5)


p3 <- ggplot(comb_table_HP1a, aes(x=protein, y=intensity_mean))

p3 <- p3 +  geom_jitter(aes(colour=protein),shape=16, alpha = 0.6, position=position_jitter(0.3)) +
  ylim(0, 100) +
  geom_boxplot(outlier.shape=NA, fill =NA) +
  scale_colour_manual(values = c("#A9A9A9", "#e33232")) +
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
# t = -9.3672, df = 163.82, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -21.14854 -13.78477
# sample estimates:
#   mean of x mean of y
# 4.443613 21.910265



################################################################################
# Plot figure

p2 + p1 + p3

filename <- "Figure_2H.pdf"
ggsave(filename,
       width = 6,
       height = 6,
       dpi = 300,
       useDingbats = FALSE)


mean(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM","intensity_mean"])
# [1] 6.701919
mean(comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 19.17108

19.17108 / 6.701919
# 2.860536



nrow(comb_table_NPM1[comb_table_NPM1$protein=="01_WT_0p5_uM",])
#[1] 180
nrow(comb_table_NPM1[comb_table_NPM1$protein=="02_MUT_0p5_uM",])
#[1] 205


mean(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM","intensity_mean"])
#[1] 21.62651
mean(comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 26.8971

26.8971 / 21.62651
# [1] 1.24371



nrow(comb_table_MED1[comb_table_MED1$protein=="01_WT_0p5_uM",])
#[1] 97
nrow(comb_table_MED1[comb_table_MED1$protein=="02_MUT_0p5_uM",])
#[1] 193





mean(comb_table_HP1a[comb_table_HP1a$protein=="01_WT_0p5_uM","intensity_mean"])
#[1] 4.443613
mean(comb_table_HP1a[comb_table_HP1a$protein=="02_MUT_0p5_uM","intensity_mean"])
# [1] 21.91026

21.91026 / 4.443613
# 4.930731




nrow(comb_table_HP1a[comb_table_HP1a$protein=="01_WT_0p5_uM",])
#[1] 143
nrow(comb_table_HP1a[comb_table_HP1a$protein=="02_MUT_0p5_uM",])
#[1] 163
