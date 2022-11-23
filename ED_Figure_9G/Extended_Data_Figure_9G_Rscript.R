################################################################################
## This script contains the code to produce Extended Data Figure 9G
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

#WT
FRAP_table_3 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP3.txt")
FRAP_table_4 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP4.txt")
FRAP_table_5 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP5.txt")
FRAP_table_6 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP6.txt")
FRAP_table_7 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP7.txt")
FRAP_table_8 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP8.txt")
FRAP_table_9 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP9.txt")
FRAP_table_10 <- read.delim(file="DVL1_FRAP_data/DVL1_WT_FRAP10.txt")

FRAP_table_3 <- normalize(FRAP_table_3)
FRAP_table_4 <- normalize(FRAP_table_4)
FRAP_table_5 <- normalize(FRAP_table_5)
FRAP_table_6 <- normalize(FRAP_table_6)
FRAP_table_7 <- normalize(FRAP_table_7)
FRAP_table_8 <- normalize(FRAP_table_8)
FRAP_table_9 <- normalize(FRAP_table_9)
FRAP_table_10 <- normalize(FRAP_table_10)


DVL1_WT <- as.data.frame(cbind(FRAP_table_3$Intensity.Region.2_NORM,
                           FRAP_table_4$Intensity.Region.2_NORM,
                           FRAP_table_5$Intensity.Region.2_NORM,
                           FRAP_table_6$Intensity.Region.2_NORM,
                           FRAP_table_7$Intensity.Region.2_NORM,
                           FRAP_table_8$Intensity.Region.2_NORM,
                           FRAP_table_9$Intensity.Region.2_NORM,
                           FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(DVL1_WT)

DVL1_WT$mean <- rowMeans(DVL1_WT[,c(1:n_columns)])
DVL1_WT$stdev <-rowSds(as.matrix(DVL1_WT[,c(1:n_columns)]))
DVL1_WT$time <- c(1:nrow(DVL1_WT))





# MUT 
FRAP_table_1 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP1.txt")
FRAP_table_2 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP2.txt")
FRAP_table_3 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP3.txt")
FRAP_table_4 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP4.txt")
FRAP_table_5 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP5.txt")
FRAP_table_6 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP6.txt")
FRAP_table_7 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP7.txt")
FRAP_table_8 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP8.txt")
FRAP_table_9 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP9.txt")
FRAP_table_10 <- read.delim(file="DVL1_FRAP_data/DVL1_MUT_FRAP10.txt")

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


DVL1_MUT <- as.data.frame(cbind(FRAP_table_1$Intensity.Region.2_NORM, 
                                     FRAP_table_2$Intensity.Region.2_NORM, 
                                     FRAP_table_3$Intensity.Region.2_NORM,
                                     FRAP_table_4$Intensity.Region.2_NORM,
                                     FRAP_table_5$Intensity.Region.2_NORM,
                                     FRAP_table_6$Intensity.Region.2_NORM,
                                     FRAP_table_7$Intensity.Region.2_NORM,
                                     FRAP_table_8$Intensity.Region.2_NORM,
                                FRAP_table_9$Intensity.Region.2_NORM,
                                FRAP_table_10$Intensity.Region.2_NORM))

n_columns <- ncol(DVL1_MUT)


DVL1_MUT$mean <- rowMeans(DVL1_MUT[,c(1:n_columns)])
DVL1_MUT$stdev <-rowSds(as.matrix(DVL1_MUT[,c(1:n_columns)]))
DVL1_MUT$time <- c(1:nrow(DVL1_MUT))

# Combined plot
DVL1_WT$protein <- "01_DVL1_WT"
DVL1_MUT$protein <- "02_DVL1_MUT"

#Select columns so that merging works
DVL1_WT_select_columns <- DVL1_WT[,c("mean","time","stdev","protein")]
DVL1_MUT_select_columns <- DVL1_MUT[,c("mean","time","stdev","protein")]

combined_table <- rbind(DVL1_WT_select_columns, DVL1_MUT_select_columns)


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

file <- "Extended_Data_Figure_9G.pdf"
ggsave(file, width = 15 , height = 10, units = "cm", useDingbats =F)

