################################################################################
## This script contains the code to produce Extended Data Figure 4J
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Shaon Basu
################################################################################

library(patchwork)

require(ggplot2)
require(RColorBrewer)
require(see)


txtfiles<- list.files(path="data", pattern="droplet.csv", full.names=T)
test<-read.csv(txtfiles[1], sep=";")
test$name <- txtfiles[1]

for(i in 2:length(txtfiles))
{
  test2 <-read.csv(txtfiles[i], sep=";")
  test2$name <- txtfiles[i]
  test <- rbind(test,test2)
}
test <- test[is.na(test$ID..ID..I)==0,]
test <- test[!duplicated(test),]
write.csv(test,"MergedTable")
test2 <- read.csv2("MergedTable", header=T, sep=",", dec=",")


test <- test2

#----------------------------------------------------------------------------------------------'

partitionGR <- data.frame(test$ImageDocumentName..Image.Name, test$IntensityMean_Ch1..Intensity.Mean.Value.of.channel..Ch1...R, test$IntensityMean_Ch2..Intensity.Mean.Value.of.channel..Ch2...R,
                          test$Area..Area..R,test$name)

condensates <- partitionGR

for(i in 1: nrow(condensates))
{
  if(grepl("_WT_", condensates[i,1])==TRUE)
  {condensates[i,6]="WT"}
  else if(grepl("_MUT_", condensates[i,1])==TRUE)
  {condensates[i,6]="MUT"}
  else
  {condensates[i,6]="NA"}
}

for(i in 1: nrow(condensates))
{
  if(grepl("MED1", condensates[i,1])==TRUE)
  {condensates[i,7]="MED1"}
  else if(grepl("NPM1", condensates[i,1])==TRUE)
  {condensates[i,7]="NPM1"}
  else if(grepl("HP1a", condensates[i,1])==TRUE)
  {condensates[i,7]="HP1a"}
}

# NOTE: Error in file naming, file name "2p5uM" is incorrect in metadata. C was 0.5 uM, this is why 2p5uM is considered as 0.5 uM.

for(i in 1: nrow(condensates))
{
  if(grepl("_2p5uM_", condensates[i,1])==TRUE)
  {condensates[i,8]=0.5}
  else if(grepl("_5uM_", condensates[i,1])==TRUE)
  {condensates[i,8]=5}
  else if(grepl("_0p5uM_", condensates[i,1])==TRUE)
  {condensates[i,8]=0.5}

}

for(i in 1: nrow(condensates))
{
  if(grepl("peptide", condensates[i,5])==TRUE)
  {condensates[i,9]='synthetic'}
  else
  {condensates[i,9]='NA'}

}


colnames(condensates) <- c('name','red','green','area','source','mix','scaffold','concentration','protein')


# Plot 0.5 uM samples
target <- condensates
target$mix <- factor(target$mix, levels = c('MUT','WT'))
target$scaffold <- factor(target$scaffold, levels = c('NPM1','MED1','HP1a'))
target <- subset(target, concentration == 0.5)

target2 <- subset(target,  mix == 'WT' | mix == 'MUT')


DFP <- ggplot(target2, aes(x=green, y=red, color = mix)) + theme_classic() + ylab('5FAM signal [HMGB1]') + xlab(paste('mCherry signal')) +
  ggtitle(paste('5FAM HMGB1 IDR, 0.5 uM')) +  scale_x_log10() + ylim(0,30) +
  scale_color_manual(values = c('red','grey40')) + geom_point(alpha = 0.6, size = target$area/5000000) +
  theme(panel.background=element_rect(fill="white", color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

DFP_plot1 <- DFP + facet_wrap(~scaffold, scales  = 'free') + theme(strip.background =element_rect(fill="white"))



# Plot 5 uM samples
target <- condensates
target$mix <- factor(target$mix, levels = c('MUT','WT'))
target$scaffold <- factor(target$scaffold, levels = c('NPM1','MED1','HP1a'))
target <- subset(target, concentration == 5)

target2 <- subset(target,  mix == 'WT' | mix == 'MUT')

DFP <- ggplot(target2, aes(x=green, y=red, color = mix)) + theme_classic() + ylab('5FAM signal [HMGB1]') + xlab(paste('mCherry signal')) +
  ggtitle(paste('5FAM HMGB1 IDR, 5 uM')) +  scale_x_log10() + ylim(0,90) +
  scale_color_manual(values = c('red','grey40')) + geom_point(alpha = 0.6, size = target$area/5000000) +
  theme(panel.background=element_rect(fill="white", color = 'black'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

DFP_plot2 <- DFP + facet_wrap(~scaffold, scales  = 'free') + theme(strip.background =element_rect(fill="white"))


# Print combined plot

DFP_plot1 / DFP_plot2

ggsave("Extended_Data_Figure_4J.pdf", width = 8, height = 10)
