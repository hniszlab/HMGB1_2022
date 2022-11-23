################################################################################
## This script contains the code to produce Figure 3D
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(matrixStats)
library(ggplot2)


normalize <- function(x) {
        x$Intensity.Region.2_NORM <- x$Intensity.Region.2/x$Intensity.Region.2[1]
        return(x)
}

################################################################################
# Read WT data

FRAP_table_1 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_1_gain_lower_final.txt")
FRAP_table_2 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_2_gain_lower.txt")
FRAP_table_3 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_3_gain_lower.txt")
FRAP_table_4 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_4_gain_lower.txt")
FRAP_table_5 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_5_gain_lower.txt")
FRAP_table_6 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_6_gain_lower.txt")
FRAP_table_7 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_7_gain_lower.txt")
FRAP_table_8 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_8_gain_lower.txt")
FRAP_table_9 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_9_gain_lower.txt")
FRAP_table_10 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_WT_FL_FIB1_FRAP_EGFP_10_gain_lower.txt")

# Using normalize function defined above

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


WT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                           FRAP_table_2$Intensity.Region.2_NORM,
                           FRAP_table_3$Intensity.Region.2_NORM,
                           FRAP_table_4$Intensity.Region.2_NORM,
                           FRAP_table_5$Intensity.Region.2_NORM,
                           FRAP_table_6$Intensity.Region.2_NORM,
                           FRAP_table_7$Intensity.Region.2_NORM,
                           FRAP_table_8$Intensity.Region.2_NORM,
                           FRAP_table_9$Intensity.Region.2_NORM,
                           FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(WT)

WT$mean <- rowMeans(WT[,c(1:n_columns)])
WT$stdev <-rowSds(as.matrix(WT[,c(1:n_columns)]))
WT$time <- c(1:nrow(WT))


p <- ggplot(WT, aes(x= time, y= mean)) +
        geom_line() + geom_point() +
        geom_errorbar(aes(ymin = mean-stdev, ymax= mean + stdev),
                      width=.2, position=position_dodge(0.05)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line()) +
        ylim(0,1.2)

p + labs(title="HMGB1 WT- EGFP FRAP", y="Signal intensity", x = "Time (s)")


################################################################################
# Read MUT data

FRAP_table_1 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_8.txt")
FRAP_table_9 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_9.txt")
FRAP_table_10 <- read.delim(file="FRAP_data_files/20210819_U2OS_HMGB1_MUT_FL_FIB1_FRAP_EGFP_gain_lower_10.txt")


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

MUT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                           FRAP_table_2$Intensity.Region.2_NORM,
                           FRAP_table_3$Intensity.Region.2_NORM,
                           FRAP_table_4$Intensity.Region.2_NORM,
                           FRAP_table_5$Intensity.Region.2_NORM,
                           FRAP_table_6$Intensity.Region.2_NORM,
                           FRAP_table_7$Intensity.Region.2_NORM,
                           FRAP_table_8$Intensity.Region.2_NORM,
                           FRAP_table_9$Intensity.Region.2_NORM,
                           FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(MUT)

MUT$mean <- rowMeans(MUT[,c(1:n_columns)])
MUT$stdev <-rowSds(as.matrix(MUT[,c(1:n_columns)]))
MUT$time <- c(1:nrow(MUT))



################################################################################
# Combined plot

WT$protein <- "WT_FL"
MUT$protein <- "MUT_FL"

WT$time <- WT$time*2
MUT$time <- MUT$time*2

WT$time <- WT$time-2
MUT$time <- MUT$time-2

combined_table <- rbind(WT, MUT)


p <- ggplot(combined_table, aes(x= time, y= mean, colour = factor(protein))) +
        geom_line() + geom_point() +
        geom_pointrange(aes(ymin = mean-stdev, ymax= mean + stdev), size = 0.2) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line()) +
        ylab("Relative signal intensity") +
        xlab("Time (s)") +
        ylim(0,1.2) +
        xlim(0,40)

p + scale_color_manual(values=c("red", "black"))

file <- "Figure_3D.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)
