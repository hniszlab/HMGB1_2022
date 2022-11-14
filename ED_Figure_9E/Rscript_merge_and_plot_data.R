# 2022.09.20

# Merge data

library(ggplot2)
library(stringr)
library(dplyr)


nucleolar_enrich_table1 <- read.delim(file="set1_print_nucleolar_enrichment_quant.txt", sep = " ")
nucleolar_enrich_table2 <- read.delim(file="set2_print_nucleolar_enrichment_quant.txt", sep = " ")
nucleolar_enrich_table3 <- read.delim(file="set3_print_nucleolar_enrichment_quant.txt", sep = " ")

nucleolar_enrich_table1 <- nucleolar_enrich_table1[,c(3,2,1)]
nucleolar_enrich_table2 <- nucleolar_enrich_table2[,c(3,2,1)]
nucleolar_enrich_table3 <- nucleolar_enrich_table3[,c(3,2,1)]

unique(nucleolar_enrich_table2$sample)


# Remove from set2 those that are part of set 3:
nucleolar_enrich_table2 <- nucleolar_enrich_table2[-grep("MEN1", nucleolar_enrich_table2$sample), ]
nucleolar_enrich_table2 <- nucleolar_enrich_table2[-grep("SQSTM1", nucleolar_enrich_table2$sample), ]
nucleolar_enrich_table2 <- nucleolar_enrich_table2[-grep("FOXL2", nucleolar_enrich_table2$sample), ]


# Rename samples (to keep intended order)
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "01_PHOX2B_WT", "13_PHOX2B_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "02_PHOX2B_MUT", "14_PHOX2B_MUT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "03_SOX2_WT", "15_SOX2_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "04_SOX2_MUT", "16_SOX2_MUT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "05_CALR_WT", "17_CALR_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "06_CALR_MUT", "18_CALR_MUT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "07_SQSTM1_WT", "19_SQSTM1_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "08_SQSTM1_MUT", "20_SQSTM1_MUT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "09_MEN1_WT", "21_MEN1_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "10_MEN1_MUT", "22_MEN1_MUT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "11_FOXL2_WT", "23_FOXL2_WT")
#nucleolar_enrich_table2$sample <- str_replace(nucleolar_enrich_table2$sample, "12_FOXL2_MUT", "24_FOXL2_MUT")

df_to_plot <- rbind(nucleolar_enrich_table1, nucleolar_enrich_table2, nucleolar_enrich_table3)


#reorder based on mutant nucleolar enrichment
df_to_plot_reorder <- df_to_plot

# Get average nucleolar enrichment
df_to_plot_reorder <- df_to_plot_reorder %>% group_by(sample) %>% mutate(group_mean = mean(log2_nucleolar_enrichment))

# need to arrange by MUT nucleolar enrichment
df_to_plot_reorder <- as.data.frame(df_to_plot_reorder)
df_to_plot_reorder$sample <- as.vector(df_to_plot_reorder$sample)

proteins <- c()
for (i in 1:nrow(df_to_plot_reorder)) {
  proteins <- append(proteins, unlist(strsplit(df_to_plot_reorder[i,1], split = "_"))[2]) 
}

df_to_plot_reorder$protein <- proteins

temp_df <- data.frame(sample = character(0), nucleolar_enrichment = double(0), 
                      log2_nucleolar_enrichment = double(0), group_mean = double(0), 
                      protein=character(0), mut_mean=double(0))

# set mut mean for each group
for (i in unique(df_to_plot_reorder$protein)) {
  
  df_subset <- df_to_plot_reorder[df_to_plot_reorder$protein==i,]
  # get mutant
  mut_df <- df_subset[grep("MUT", df_subset$sample), ]
  df_subset$mut_mean <- mut_df$group_mean[1]
  temp_df <- rbind(temp_df, df_subset)
}

df_to_plot_reordered <- temp_df

df_to_plot_reordered$sample <- as.factor(df_to_plot_reordered$sample)
df_to_plot_reordered$sample <- reorder(df_to_plot_reordered$sample, -df_to_plot_reordered$mut_mean)



# Plot data as boxplot + jitter

p <- ggplot(df_to_plot_reordered, aes(x=sample, y=log2_nucleolar_enrichment))

p +  geom_jitter(aes(colour=sample),shape=16, alpha = 0.8, position=position_jitter(0.3)) +
  #ylim(0, 30) + 
  geom_boxplot(outlier.shape=NA, fill =NA) +
  #stat_summary(fun=mean, geom="point", shape=18,
  #                         size=3, color="red") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  scale_colour_manual(values = rep(c("#A9A9A9","#D22B2B"),12)) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 90))
NULL

file <- "Extended_Data_Figure_9E.pdf"
ggsave(file, width = 24 , height = 15, units = "cm")





# Check n of measurements

nrow(df_to_plot[grep("HMGB3_WT", df_to_plot$sample), ])
# 24
nrow(df_to_plot[grep("HMGB3_MUT", df_to_plot$sample), ])
# 20

nrow(df_to_plot[grep("FOXC1_WT", df_to_plot$sample), ])
# 21
nrow(df_to_plot[grep("FOXC1_MUT", df_to_plot$sample), ])
# 28

nrow(df_to_plot[grep("PHOX2B_WT", df_to_plot$sample), ])
# 41
nrow(df_to_plot[grep("PHOX2B_MUT", df_to_plot$sample), ])
# 33

nrow(df_to_plot[grep("RAX_WT", df_to_plot$sample), ])
# 31
nrow(df_to_plot[grep("RAX_MUT", df_to_plot$sample), ])
# 35

nrow(df_to_plot[grep("MYOD_WT", df_to_plot$sample), ])
# 33
nrow(df_to_plot[grep("MYOD_MUT", df_to_plot$sample), ])
# 30

nrow(df_to_plot[grep("MEN1_WT", df_to_plot$sample), ])
# 63
nrow(df_to_plot[grep("MEN1_MUT", df_to_plot$sample), ])
# 37

nrow(df_to_plot[grep("FOXF1_WT", df_to_plot$sample), ])
# 17
nrow(df_to_plot[grep("FOXF1_MUT", df_to_plot$sample), ])
# 17

nrow(df_to_plot[grep("SOX2_WT", df_to_plot$sample), ])
# 72
nrow(df_to_plot[grep("SOX2_MUT", df_to_plot$sample), ])
# 57

nrow(df_to_plot[grep("CALR_WT", df_to_plot$sample), ])
# 83
nrow(df_to_plot[grep("CALR_MUT", df_to_plot$sample), ])
# 36

nrow(df_to_plot[grep("SQSTM1_WT", df_to_plot$sample), ])
# 47
nrow(df_to_plot[grep("SQSTM1_MUT", df_to_plot$sample), ])
# 53

nrow(df_to_plot[grep("RUNX1_WT", df_to_plot$sample), ])
# 28
nrow(df_to_plot[grep("RUNX1_MUT", df_to_plot$sample), ])
# 31

nrow(df_to_plot[grep("FOXL2_WT", df_to_plot$sample), ])
# 30
nrow(df_to_plot[grep("FOXL2_MUT", df_to_plot$sample), ])
# 69



# Print values:
df_to_plot <- rbind(nucleolar_enrich_table1, nucleolar_enrich_table2, nucleolar_enrich_table3)

write.table(df_to_plot, file ="set1_set2_set3_print_nucleolar_enrichment_quant.txt", row.names = F)

# Print table for correlation plot, exclude SQSTM1 WT:
print_table_of_values2 <- df_to_plot[-grep("SQSTM1_WT", df_to_plot$sample), ]
write.table(print_table_of_values2, file ="set1_set2_set3_print_nucleolar_enrichment_quant_noSQSTM1_WT.txt", row.names = F)


