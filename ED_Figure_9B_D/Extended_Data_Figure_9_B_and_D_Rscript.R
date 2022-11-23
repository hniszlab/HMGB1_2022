################################################################################
## This script contains the code to produce Extended Data Figure 9B and 9D
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Henri Niskanen
################################################################################

library(ggplot2)

################################################################################
# HMGB1

charge_result <- read.delim(file="charge_data/Charge_HMGB1wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_HMGB1mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"


df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "HMGB1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_HMGB1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# HMGB3

charge_result <- read.delim(file="charge_data/Charge_HMGB3_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_HMGB3_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"


df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "HMGB3 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_HMGB3.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# FOXC1

charge_result <- read.delim(file="charge_data/Charge_FOXC1wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_FOXC1mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXC1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_FOXC1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# FOXF1

charge_result <- read.delim(file="charge_data/Charge_FOXF1wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_FOXF1mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"


df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXF1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_FOXF1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# MYOD1

charge_result <- read.delim(file="charge_data/Charge_MYOD1wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_MYOD1mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "MYOD1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_MYOD1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# RAX

charge_result <- read.delim(file="charge_data/Charge_RAXwt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_RAXmut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "RAX Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_RAX.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)


################################################################################
# RUNX1

charge_result <- read.delim(file="charge_data/Charge_RUNX1wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_RUNX1mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "RUNX1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_B_RUNX1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# CALR

charge_result <- read.delim(file="charge_data/Charge_CALR_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_CALR_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "CALR Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_CALR.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# FOXL2

charge_result <- read.delim(file="charge_data/Charge_FOXL2_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_FOXL2_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "FOXL2 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_FOXL2.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# PHOX2B

charge_result <- read.delim(file="charge_data/Charge_PHOX2B_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_PHOX2B_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "PHOX2B Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_PHOX2B.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# SOX2

charge_result <- read.delim(file="charge_data/Charge_SOX2_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_SOX2_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"


df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "SOX2 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_SOX2.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)

################################################################################
# SQSTM1

charge_result <- read.delim(file="charge_data/Charge_SQSTM1_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_SQSTM1_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "SQSTM1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_SQSTM1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# MEN1

charge_result <- read.delim(file="charge_data/Charge_MEN1_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_MEN1_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "MEN1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_MEN1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)



################################################################################
# DVL1

charge_result <- read.delim(file="charge_data/Charge_DVL1_wt.txt")
colnames(charge_result) <- c("Position","Residue","Charge")
charge_result$Position <- c(1:nrow(charge_result))
charge_result$type <- "WT"

charge_result2 <- read.delim(file="charge_data/Charge_DVL1_mut.txt")
colnames(charge_result2) <- c("Position","Residue","Charge")
charge_result2$Position <- c(1:nrow(charge_result2))
charge_result2$type <- "MUT"

df <- rbind(charge_result, charge_result2)

ggplot(data=df, aes(x=Position, y=Charge, group = type)) +
  geom_line(aes(color = type, linetype=type), size=1.5) +
  labs(y ="VSL2", x = "AA residue", title = "DVL1 Charge plot") +
  scale_color_manual(values=c('red','#58595b')) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  theme_classic() +
  geom_hline(yintercept= 0, linetype="dotted", color = "black", size=1) +
  ylim(-1,1)

file <- "Extended_Data_Figure_9_D_DVL1.pdf"

ggsave(file, width = 15 , height = 13, units = "cm", useDingbats =F)
