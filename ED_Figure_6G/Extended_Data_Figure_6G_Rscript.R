################################################################################
## This script contains the code to produce Extended Data Figure 6G
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(matrixStats)
library(ggplot2)

FRAP_table_1 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP8.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP9.txt")
FRAP_table_10 <- read.delim(file="FRAP_data/20220408_PBmix_WT_DOX_24h_FRAP10.txt")

normalize <- function(x) {
  x$Intensity.Region.2_NORM <- x$Intensity.Region.2/x$Intensity.Region.2[1]
  return(x)
}

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

B1_WT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                           FRAP_table_2$Intensity.Region.2_NORM,
                           FRAP_table_3$Intensity.Region.2_NORM,
                           FRAP_table_4$Intensity.Region.2_NORM,
                           FRAP_table_5$Intensity.Region.2_NORM,
                           FRAP_table_6$Intensity.Region.2_NORM,
                           FRAP_table_7$Intensity.Region.2_NORM,
                           FRAP_table_8$Intensity.Region.2_NORM,
                           FRAP_table_9$Intensity.Region.2_NORM,
                           FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(B1_WT)

B1_WT$mean <- rowMeans(B1_WT[,c(1:n_columns)])
B1_WT$stdev <-rowSds(as.matrix(B1_WT[,c(1:n_columns)]))
B1_WT$time <- c(1:nrow(B1_WT))

# MUT
FRAP_table_1 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP8.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP9.txt")
FRAP_table_10 <- read.delim(file="FRAP_data/20220408_PBmix_MUT_DOX_24h_FRAP10.txt")

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

B1_MUT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                                  FRAP_table_2$Intensity.Region.2_NORM,
                                  FRAP_table_3$Intensity.Region.2_NORM,
                                  FRAP_table_4$Intensity.Region.2_NORM,
                                  FRAP_table_5$Intensity.Region.2_NORM,
                                  FRAP_table_6$Intensity.Region.2_NORM,
                                  FRAP_table_7$Intensity.Region.2_NORM,
                                  FRAP_table_8$Intensity.Region.2_NORM,
                                  FRAP_table_9$Intensity.Region.2_NORM,
                                  FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT)

B1_MUT$mean <- rowMeans(B1_MUT[,c(1:n_columns)])
B1_MUT$stdev <-rowSds(as.matrix(B1_MUT[,c(1:n_columns)]))
B1_MUT$time <- c(1:nrow(B1_MUT))



# MUT delTail
FRAP_table_1 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP1.txt")
FRAP_table_2 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP2.txt")
FRAP_table_3 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP3.txt")
FRAP_table_4 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP4.txt")
FRAP_table_5 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP5.txt")
FRAP_table_6 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP6.txt")
FRAP_table_7 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP7.txt")
FRAP_table_8 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP8.txt")
FRAP_table_9 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP9.txt")
FRAP_table_10 <- read.delim(file="FRAP_data/20220408_PBmix_MUTdeltail_DOX_24h_FRAP10.txt")

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

B1_MUT_deltail <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM,
                              FRAP_table_2$Intensity.Region.2_NORM,
                              FRAP_table_3$Intensity.Region.2_NORM,
                              FRAP_table_4$Intensity.Region.2_NORM,
                              FRAP_table_5$Intensity.Region.2_NORM,
                              FRAP_table_6$Intensity.Region.2_NORM,
                              FRAP_table_7$Intensity.Region.2_NORM,
                              FRAP_table_8$Intensity.Region.2_NORM,
                              FRAP_table_9$Intensity.Region.2_NORM,
                              FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(B1_MUT_deltail)

B1_MUT_deltail$mean <- rowMeans(B1_MUT_deltail[,c(1:n_columns)])
B1_MUT_deltail$stdev <-rowSds(as.matrix(B1_MUT_deltail[,c(1:n_columns)]))
B1_MUT_deltail$time <- c(1:nrow(B1_MUT_deltail))


# Combine tables and plot data:

B1_WT$protein <- "01_HMGB1_WT"
B1_MUT$protein <- "02_HMGB1_MUT"
B1_MUT_deltail$protein <- "03_HMGB1_MUT_Patchless"

B1_WT_select_columns <- B1_WT[,c("mean","time","stdev","protein")]
B1_MUT_select_columns <- B1_MUT[,c("mean","time","stdev","protein")]
B1_MUT_deltail_select_columns <- B1_MUT_deltail[,c("mean","time","stdev","protein")]

combined_table <- rbind(B1_WT_select_columns,
                        B1_MUT_select_columns,
                        B1_MUT_deltail_select_columns)

p <- ggplot(combined_table, aes(x= time, y= mean, colour = factor(protein))) +
  geom_line() + geom_point() +
  geom_pointrange(aes(ymin = mean-stdev, ymax= mean + stdev), size = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  ylim(0,1.2) +
  xlim(0,40)

p + scale_color_manual(values=c("black", "red", "purple")) +
  annotate("text", x = 30, y = 0.20, colour = "black", label = "n = 10") +
  annotate("text", x = 30, y = 0.15, colour = "red", label = "n = 10") +
  annotate("text", x = 30, y = 0.10, colour = "purple", label = "n = 10")

file <- "Extended_Data_Figure_6G.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)
