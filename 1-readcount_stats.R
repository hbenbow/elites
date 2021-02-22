library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)
allowWGCNAThreads()

setwd("~/Documents/Kat/Elites/")
# where are kalliso files?
load("Elite_Panel.RData")
# Import kallisto files with txi 
# ==================================================================================
colnames(colData)<-c("Code", "Rep", "Timepoint", "Genotype", "Treatment", "Sample", "Factor")
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)


# check that order of samples in metadata and txi object are the same

# ==================================================================================
# read count stats chunk starts here

expressed_genes<-txi.kallisto.tsv$abundance
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
# expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),] #only use if you dont want LC genes
expressed_genes<-expressed_genes[,c(49, 1:48)]
expressed_genes_long<-expressed_genes %>% gather(Code, TPM, 2:49)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Code")
sub<-all_wheat_genes[,c(9, 2, 3, 4)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>=0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>=0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>=0.5, 1,0)
rep_wise$Sum<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3

expressed<-rep_wise[(rep_wise$Sum >=2),]

for(i in unique(expressed$Factor)){
  data<-expressed[(expressed$Factor==i),]
  factor<-paste(i)
  write.csv(data, file=paste("~/Documents/Kat/Elites/data/", factor, ".csv", sep=""))
  assign(factor, data)
}
{
  R1<-rbind(C1R, T1R)
  R2<-rbind(C2R, T2R)
  R3<-rbind(C3R, T3R)
  R4<-rbind(C4R, T4R)
  R1<-R1[!(duplicated(R1$GeneID)),]
  R2<-R2[!(duplicated(R2$GeneID)),]
  R3<-R3[!(duplicated(R3$GeneID)),]
  R4<-R4[!(duplicated(R4$GeneID)),]
  
  S1<-rbind(C1S, T1S)
  S2<-rbind(C2S, T2S)
  S3<-rbind(C3S, T3S)
  S4<-rbind(C4S, T4S)
  S1<-S1[!(duplicated(S1$GeneID)),]
  S2<-S2[!(duplicated(S2$GeneID)),]
  S3<-S3[!(duplicated(S3$GeneID)),]
  S4<-S4[!(duplicated(S4$GeneID)),]
  
  

  R1$Comparison<-"R1"
  R2$Comparison<-"R2"
  R3$Comparison<-"R3"
  R4$Comparison<-"R4"
  S1$Comparison<-"S1"
  S2$Comparison<-"S2"
  S3$Comparison<-"S3"
  S4$Comparison<-"S4"
  
  R1<-R1[,c(2, 10)]
  R2<-R2[,c(2, 10)]
  R3<-R3[,c(2, 10)]
  R4<-R4[,c(2, 10)]
  S1<-S1[,c(2, 10)]
  S2<-S2[,c(2, 10)]
  S3<-S3[,c(2, 10)]
  S4<-S4[,c(2, 10)]
  }

all_filtered_lists<-rbind(R1,
                          R2,
                          R3,
                          R4,
                          S1,
                          S2,
                          S3,
                          S4)

kable(table(all_filtered_lists$Comparison))

write.csv(all_filtered_lists, file="~/Documents/Kat/Elites/data/all_filtered_lists.csv", row.names = F)
# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Genotype", "Rep", "Timepoint"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Timepoint), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black"))

