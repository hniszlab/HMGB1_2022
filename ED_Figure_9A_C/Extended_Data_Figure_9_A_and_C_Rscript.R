################################################################################
## This script contains the code to produce Extended Data Figure 9A and 9C
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)

################################################################################
# HMGB1

PONDR_result <- read.delim(file="PONDR_values/HMGB1_WT_FL.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/HMGB1_MUT_FL.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VLXT, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VLXT", x = "AA residue", title = "HMGB1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_HMGB1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# HMGB3

PONDR_result <- read.delim(file="PONDR_values/HMGB3_wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/HMGB3_mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VLXT, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VLXT", x = "AA residue", title = "HMGB3 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_HMGB3.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)


################################################################################
# FOXC1

PONDR_result <- read.delim(file="PONDR_values/FOXC1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/FOXC1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXC1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_FOXC1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# FOXF1

PONDR_result <- read.delim(file="PONDR_values/FOXF1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/FOXF1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXF1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_FOXF1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)





################################################################################
# MYOD1

PONDR_result <- read.delim(file="PONDR_values/MYOD1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/MYOD1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "MYOD1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_MYOD1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)


################################################################################
# RAX

PONDR_result <- read.delim(file="PONDR_values/RAXwt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/RAXmut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "RAX PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9A_RAX.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)




################################################################################
# RUNX1

PONDR_result <- read.delim(file="PONDR_values/RUNX1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/RUNX1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "RUNX1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_RUNX1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# CALR

PONDR_result <- read.delim(file="PONDR_values/PONDR_CALR_wt.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/PONDR_CALR_mut.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "CALR PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_CALR.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)




################################################################################
# FOXL2

PONDR_result <- read.delim(file="PONDR_values/PONDR_FOXL2_wt.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/PONDR_FOXL2_mut.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXL2 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_FOXL2.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)





################################################################################
# PHOX2B

PONDR_result <- read.delim(file="PONDR_values/PONDR_PHOX2B_wt.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/PONDR_PHOX2B_mut.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "PHOX2B PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_PHOX2B.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)




################################################################################
# SOX2

PONDR_result <- read.delim(file="PONDR_values/PONDR_SOX2_wt.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/PONDR_SOX2_mut.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "SOX2 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_SOX2.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# SQSTM1

PONDR_result <- read.delim(file="PONDR_values/PONDR_SQSTM1wt.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/PONDR_SQSTM1_mut.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "SQSTM1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_SQSTM1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# MEN1

PONDR_result <- read.delim(file="PONDR_values/MEN1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/MEN1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "MEN1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_MEN1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)




################################################################################
# DVL1

PONDR_result <- read.delim(file="PONDR_values/DVL1wt_PONDR.txt")
colnames(PONDR_result) <- c("Num","Res","VLXT","VSL2","type")
PONDR_result$type <- "WT"

PONDR_result2 <- read.delim(file="PONDR_values/DVL1mut_PONDR.txt")
colnames(PONDR_result2) <- c("Num","Res","VLXT","VSL2", "type")
PONDR_result2$type <- "MUT"

df <- rbind(PONDR_result, PONDR_result2)

ggplot(data=df, aes(x=Num, y=VSL2, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "DVL1 PONDR PLOT") +
  scale_color_manual(values=c('red','black')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0.5, linetype="dotted", color = "black", size=1) +
  ylim(0,1)

file <- "Extended_Data_Figure_9C_DVL1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)
