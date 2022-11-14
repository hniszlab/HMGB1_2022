################################################################################
## This script contains the code to produce Figure 3G
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(matrixStats)
library(ggplot2)
#library(RColorBrewer)

normalize <- function(x) {
  x$Intensity.Region.2_NORM <- x$Intensity.Region.2/x$Intensity.Region.2[1]
  return(x)
}

FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_4.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_5.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_6.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_7.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_8.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_9.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/HMGB1_WT_FRAP_10.txt")


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

# MUT
FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_MUT_FRAP_8.txt")

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

#Select columns so that merging works
B1_WT_select_columns <- B1_WT[,c("mean","time","stdev","protein")]
B1_MUT_select_columns <- B1_MUT[,c("mean","time","stdev","protein")]

combined_table <- rbind(B1_WT_select_columns, B1_MUT_select_columns)


p <- ggplot(data = combined_table, aes(x = time, group = protein)) +
  geom_line(aes(y = mean, color = protein), size = 1) +
  geom_ribbon(aes(y = mean, ymin = mean - stdev, ymax = mean + stdev, fill = protein), alpha = .2) +
  xlab("Time (s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())



# Rdel

FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_6.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_7.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_8.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_9.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/HMGB1_delFS_FRAP_10.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)

B1_delFS<- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                   FRAP_table_2$Intensity.Region.2_NORM,
                                   FRAP_table_3$Intensity.Region.2_NORM,
                                   FRAP_table_4$Intensity.Region.2_NORM,
                                   FRAP_table_5$Intensity.Region.2_NORM,
                                   FRAP_table_6$Intensity.Region.2_NORM,
                                FRAP_table_7$Intensity.Region.2_NORM,
                               FRAP_table_8$Intensity.Region.2_NORM,
                               FRAP_table_9$Intensity.Region.2_NORM))


n_columns <- ncol(B1_delFS)

B1_delFS$mean <- rowMeans(B1_delFS[,c(1:n_columns)])
B1_delFS$stdev <-rowSds(as.matrix(B1_delFS[,c(1:n_columns)]))
B1_delFS$time <- c(1:nrow(B1_delFS))

B1_delFS$protein <- "03_HMGB1_delFS"
B1_delFS_select_columns <- B1_delFS[,c("mean","time","stdev","protein")]




##### RtoK
FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoK_FRAP_6.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)

B1_MUT_RtoK <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                 FRAP_table_2$Intensity.Region.2_NORM,
                                 FRAP_table_3$Intensity.Region.2_NORM,
                                 FRAP_table_4$Intensity.Region.2_NORM,
                                 FRAP_table_5$Intensity.Region.2_NORM,
                                 FRAP_table_6$Intensity.Region.2_NORM))


n_columns <- ncol(B1_MUT_RtoK)

B1_MUT_RtoK$mean <- rowMeans(B1_MUT_RtoK[,c(1:n_columns)])
B1_MUT_RtoK$stdev <-rowSds(as.matrix(B1_MUT_RtoK[,c(1:n_columns)]))
B1_MUT_RtoK$time <- c(1:nrow(B1_MUT_RtoK))

B1_MUT_RtoK$protein <- "04_HMGB1_MUT_RtoK"
B1_MUT_RtoK_select_columns <- B1_MUT_RtoK[,c("mean","time","stdev","protein")]



##### KtoR
FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_MUT_KtoR_FRAP_8.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)

B1_MUT_KtoR <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                 FRAP_table_2$Intensity.Region.2_NORM,
                                 FRAP_table_3$Intensity.Region.2_NORM,
                                 FRAP_table_4$Intensity.Region.2_NORM,
                                 FRAP_table_5$Intensity.Region.2_NORM,
                                 FRAP_table_6$Intensity.Region.2_NORM,
                                 FRAP_table_7$Intensity.Region.2_NORM,
                                 FRAP_table_8$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT_KtoR)
B1_MUT_KtoR$mean <- rowMeans(B1_MUT_KtoR[,c(1:n_columns)])
B1_MUT_KtoR$stdev <-rowSds(as.matrix(B1_MUT_KtoR[,c(1:n_columns)]))
B1_MUT_KtoR$time <- c(1:nrow(B1_MUT_KtoR))

B1_MUT_KtoR$protein <- "05_HMGB1_MUT_KtoR"
B1_MUT_KtoR_select_columns <- B1_MUT_KtoR[,c("mean","time","stdev","protein")]

##### R to A

FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_8.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/HMGB1_MUT_RtoA_FRAP_9.txt")

FRAP_table_1 <- normalize(FRAP_table_1)
FRAP_table_2 <- normalize(FRAP_table_2)
FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)
FRAP_table_9 <- normalize(FRAP_table_9)

B1_MUT_RtoA <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                   FRAP_table_2$Intensity.Region.2_NORM,
                                   FRAP_table_3$Intensity.Region.2_NORM,
                                   FRAP_table_4$Intensity.Region.2_NORM,
                                   FRAP_table_5$Intensity.Region.2_NORM,
                                   FRAP_table_6$Intensity.Region.2_NORM,
                                   FRAP_table_7$Intensity.Region.2_NORM,
                                   FRAP_table_8$Intensity.Region.2_NORM,
                                   FRAP_table_9$Intensity.Region.2_NORM))


n_columns <- ncol(B1_MUT_RtoA)

B1_MUT_RtoA$mean <- rowMeans(B1_MUT_RtoA[,c(1:n_columns)])
B1_MUT_RtoA$stdev <-rowSds(as.matrix(B1_MUT_RtoA[,c(1:n_columns)]))
B1_MUT_RtoA$time <- c(1:nrow(B1_MUT_RtoA))

# Combined plot

B1_MUT_RtoA$protein <- "06_HMGB1_MUT_RtoA"
B1_MUT_RtoA_select_columns <- B1_MUT_RtoA[,c("mean","time","stdev","protein")]

##### MUT del Tail
FRAP_table_1 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_7.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_8.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_9.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_10.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_11.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_12.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_13.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_14.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_15.txt")
FRAP_table_10 <- read.delim(file="FRAP_data/HMGB1_MUT_delTail_95p_FRAP_16.txt")

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

B1_MUT_delTail<- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                     FRAP_table_2$Intensity.Region.2_NORM,
                                     FRAP_table_3$Intensity.Region.2_NORM,
                                     FRAP_table_4$Intensity.Region.2_NORM,
                                     FRAP_table_5$Intensity.Region.2_NORM,
                                     FRAP_table_6$Intensity.Region.2_NORM,
                                     FRAP_table_7$Intensity.Region.2_NORM,
                                     FRAP_table_8$Intensity.Region.2_NORM,
                                     FRAP_table_9$Intensity.Region.2_NORM,
                                     FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT_delTail)

B1_MUT_delTail$mean <- rowMeans(B1_MUT_delTail[,c(1:n_columns)])
B1_MUT_delTail$stdev <-rowSds(as.matrix(B1_MUT_delTail[,c(1:n_columns)]))
B1_MUT_delTail$time <- c(1:nrow(B1_MUT_delTail))

B1_MUT_delTail$protein <- "07_HMGB1_MUT_delTail"
B1_MUT_delTail_select_columns <- B1_MUT_delTail[,c("mean","time","stdev","protein")]


# Combine all data frames
combined_table <- rbind(B1_WT_select_columns, B1_MUT_select_columns, B1_delFS_select_columns,
                        B1_MUT_RtoK_select_columns, B1_MUT_KtoR_select_columns, B1_MUT_RtoA_select_columns,
                        B1_MUT_delTail_select_columns)

# Plot data

p <- ggplot(data = combined_table, aes(x = time, group = protein)) +
  geom_line(aes(y = mean, color = protein), size = 1) +
  xlab("Time (s)") +
  ylim(0,1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  geom_ribbon(aes(y = mean, ymin = mean - stdev, ymax = mean + stdev, fill = protein), alpha = .2) +
  NULL

p

file <- "Figure_3G.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)
