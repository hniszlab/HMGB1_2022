################################################################################
## This script contains the code to produce Extended Data Figure 5C
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

# rowSds from matrixStats package:
library(matrixStats)
library(ggplot2)

FRAP_table_1 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_4.txt")
FRAP_table_4 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_5.txt")
FRAP_table_5 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_6.txt")
FRAP_table_6 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_7.txt")
FRAP_table_7 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_8.txt")
FRAP_table_8 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_9.txt")
FRAP_table_9 <- read.delim(file="FRAP_analysis_files/HMGB1_WT_FRAP_10.txt")

normalize <- function(x) {
  x$Intensity.Region.2_NORM <- x$Intensity.Region.2/x$Intensity.Region.2[1]
  return(x)
}

# WT
FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)
FRAP_table_9 <- normalize(FRAP_table_9)

B1_WT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                           FRAP_table_2$Intensity.Region.2_NORM,
                           FRAP_table_3$Intensity.Region.2_NORM,
                           FRAP_table_4$Intensity.Region.2_NORM,
                           FRAP_table_5$Intensity.Region.2_NORM,
                           FRAP_table_6$Intensity.Region.2_NORM,
                           FRAP_table_7$Intensity.Region.2_NORM,
                           FRAP_table_8$Intensity.Region.2_NORM,
                           FRAP_table_9$Intensity.Region.2_NORM))

n_columns <- ncol(B1_WT)

B1_WT$mean <- rowMeans(B1_WT[,c(1:n_columns)])
B1_WT$stdev <-rowSds(as.matrix(B1_WT[,c(1:n_columns)]))
B1_WT$time <- c(1:nrow(B1_WT))

# Mutant
FRAP_table_1 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_FRAP_8.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)

B1_MUT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                     FRAP_table_2$Intensity.Region.2_NORM,
                                     FRAP_table_3$Intensity.Region.2_NORM,
                                     FRAP_table_4$Intensity.Region.2_NORM,
                                     FRAP_table_5$Intensity.Region.2_NORM,
                                     FRAP_table_6$Intensity.Region.2_NORM,
                                     FRAP_table_7$Intensity.Region.2_NORM,
                                     FRAP_table_8$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT)

B1_MUT$mean <- rowMeans(B1_MUT[,c(1:n_columns)])
B1_MUT$stdev <-rowSds(as.matrix(B1_MUT[,c(1:n_columns)]))
B1_MUT$time <- c(1:nrow(B1_MUT))

# Combined plot
B1_WT$protein <- "01_HMGB1_WT"
B1_MUT$protein <- "02_HMGB1_MUT"

#Select columns
B1_WT_select_columns <- B1_WT[,c("mean","time","stdev","protein")]
B1_MUT_select_columns <- B1_MUT[,c("mean","time","stdev","protein")]


#### IDR1
FRAP_table_1 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_8.txt")
FRAP_table_9 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_9.txt")
FRAP_table_10 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR1_FRAP_10.txt")


FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)
FRAP_table_9 <- normalize(FRAP_table_9)
FRAP_table_10 <- normalize(FRAP_table_10)

B1_MUT_IDR1 <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                   FRAP_table_2$Intensity.Region.2_NORM,
                                   FRAP_table_3$Intensity.Region.2_NORM,
                                   FRAP_table_4$Intensity.Region.2_NORM,
                                   FRAP_table_5$Intensity.Region.2_NORM,
                                   FRAP_table_6$Intensity.Region.2_NORM,
                                   FRAP_table_7$Intensity.Region.2_NORM,
                                   FRAP_table_8$Intensity.Region.2_NORM,
                                   FRAP_table_9$Intensity.Region.2_NORM,
                                   FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT_IDR1)

B1_MUT_IDR1$mean <- rowMeans(B1_MUT_IDR1[,c(1:n_columns)])
B1_MUT_IDR1$stdev <-rowSds(as.matrix(B1_MUT_IDR1[,c(1:n_columns)]))
B1_MUT_IDR1$time <- c(1:nrow(B1_MUT_IDR1))



B1_MUT_IDR1$protein <- "03_HMGB1_MUT_IDR1"
B1_MUT_IDR1_select_columns <- B1_MUT_IDR1[,c("mean","time","stdev","protein")]


# IDR 2
FRAP_table_1 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_7.txt")
FRAP_table_7 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_8.txt")
FRAP_table_8 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_9.txt")
FRAP_table_9 <- read.delim(file="FRAP_analysis_files/HMGB1_MUT_IDR2_FRAP_10.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)
FRAP_table_9 <- normalize(FRAP_table_9)

B1_MUT_IDR2 <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                   FRAP_table_2$Intensity.Region.2_NORM,
                                   FRAP_table_3$Intensity.Region.2_NORM,
                                   FRAP_table_4$Intensity.Region.2_NORM,
                                   FRAP_table_5$Intensity.Region.2_NORM,
                                   FRAP_table_6$Intensity.Region.2_NORM,
                                   FRAP_table_7$Intensity.Region.2_NORM,
                                   FRAP_table_8$Intensity.Region.2_NORM,
                                   FRAP_table_9$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT_IDR2)

B1_MUT_IDR2$mean <- rowMeans(B1_MUT_IDR2[,c(1:n_columns)])
B1_MUT_IDR2$stdev <-rowSds(as.matrix(B1_MUT_IDR2[,c(1:n_columns)]))
B1_MUT_IDR2$time <- c(1:nrow(B1_MUT_IDR2))

# Combined plot
B1_MUT_IDR2$protein <- "04_HMGB1_MUT_IDR2"

B1_MUT_IDR2_select_columns <- B1_MUT_IDR2[,c("mean","time","stdev","protein")]

# Combined plot with WT and MUT and IDR 1 and 2
combined_table <- rbind(B1_WT_select_columns, B1_MUT_select_columns, B1_MUT_IDR1_select_columns, B1_MUT_IDR2_select_columns)

p <- ggplot(data = combined_table, aes(x = time, group = protein)) +
  geom_line(aes(y = mean, color = protein), size = 1) +
  geom_ribbon(aes(y = mean, ymin = mean - stdev, ymax = mean + stdev, fill = protein), alpha = .2) +
  ylim(0,1.2) +
  xlab("Time (s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
p

file <- "Extended_Data_Figure_5C.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)
