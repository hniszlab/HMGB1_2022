################################################################################
## This script contains the code to quantify imaging results for Figure 3I
## Mensah & Niskanen et al.
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022
## Author: Rene Buschow
################################################################################

files<-list.files(pattern = "ZOI.csv")

##### simple linear thresholds

x_cutoff1 <- 5
x_cutoff2 <- 100

###################

xpar<-4

for (n in 1:length(files))

{

  colnames<-c("Condition", "Region","Mean_Cy2_neg","Mean_Puro_neg","Cor Cy2 vs Cy5 neg","Cor p-value pos","Mean_Cy2_pos","Mean_Puro_pos","Cor Cy2 vs Cy5 pos","Cor p-value pos")
  result<-matrix(nrow=100,ncol=10)
  data_raw<-read.csv(files[n], header = FALSE, skip=2, dec=",", sep=";")


  condi<-levels(as.factor(data_raw[,13]))



  for (l in 1:length(condi))
  {

    data_raw1<-subset(data_raw,(data_raw[,13] ==  condi[l]))

    region<-levels(as.factor(data_raw1[,14]))

    for (m in 1:length(region))
    {

      data_work<-subset(data_raw1,(data_raw1[,14] ==  region[m]))

      data_gfp_neg<-subset(data_work,(data_work[,xpar] <  x_cutoff1) )
      data_gfp_pos<-subset(data_work,(data_work[,xpar] >  x_cutoff1) & (data_work[,xpar] <   x_cutoff2))

      tiff(paste(substr(files[n],1,5),"_",condi[l],"_",region[m],".tiff", sep=""),width = 500, height = 500, units="px", pointsize = 20)
      smoothScatter(data_work[,4],data_work[,5],
                    xlab="ZOI Mean GFP",
                    ylab="ZOI Mean Puro",
                    main=paste(substr(files[n],1,5),"_",condi[l],"_",region[m],sep=""), xlim=c(0,50),ylim=c(0,50),
                    nrpoints=100,
                    transformation=function(x)x/0.8,
                    bandwidth=0.4,
                    nbin=250,
                    colramp=colorRampPalette(c("white",rainbow(12))))
      dev.off()



      result[2*l+m,1]<-condi[l]
      result[2*l+m,2]<-region[m]



      if (length(data_gfp_neg[,5])>78)
      {
        spearcor<-cor.test(data_gfp_neg[,4],data_gfp_neg[,5], methode="spearman",na.action="na.omit")
      }
      else{}

      result[2*l+m,3]<-mean(data_gfp_neg[,4])
      result[2*l+m,4]<-mean(data_gfp_neg[,5])
      result[2*l+m,5]<-as.numeric(spearcor[4])
      result[2*l+m,6]<-as.numeric(spearcor[3])


      if (length(data_gfp_pos[,1])>78)
      {
        spearcor<-cor.test(data_gfp_pos[,4],data_gfp_pos[,5], methode="spearman",na.action="na.omit")
      }
      else{}

      result[2*l+m,7]<-mean(data_gfp_pos[,4])
      result[2*l+m,8]<-mean(data_gfp_pos[,5])
      result[2*l+m,9]<-as.numeric(spearcor[4])
      result[2*l+m,10]<-as.numeric(spearcor[3])



    }
  }

  write.table(result,paste(files[n],"_cor_.txt",sep=""), col.names=colnames, row.names=FALSE, sep="\t")
  rm(result)
}
