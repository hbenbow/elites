# de testing for individual datasets that were filtered before

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
library(knitr)

setwd("data/")

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files (already made so load txi object)
load("~/Documents/Kat/Elites/tximpost_all.RData")
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)

filter<-read.csv(file="~/Documents/Kat/Elites/data/all_filtered_lists.csv")
counts<-txi.kallisto.tsv$counts
tpm<-txi.kallisto.tsv$abundance

list<-list()
degs_filtered<-list()
degs_no_filter<-list()
for(set in unique(filter$Comparison)){
  comp<-filter[(filter$Comparison==set),]
  genotype<-strsplit(set, split = "")[[1]][1]
  timepoint<-strsplit(set, split = "")[[1]][2]
  metadata<-colData[(colData$Genotype==genotype),,drop=T]
  metadata<-metadata[(metadata$Timepoint==timepoint),,drop=T]
  metadata<-droplevels(metadata)
  samples<-as.character(metadata$Code)
  # genes<-as.character(comp$GeneID)
  df<-counts[,samples] #put genes back here
  case<-colnames(df)==metadata$Code
  list[[length(list)+1]]<-case
  dds <- DESeqDataSetFromMatrix(round(df), metadata, ~ Treatment)
  dds<-DESeq(dds)
  r<-as.data.frame(results(dds, contrast=c("Treatment", "T", "C") ,format='DataFrame', tidy=TRUE))
  r$comparison<-paste(set)
  degs_filtered[[length(degs_filtered)+1]]<-r
  
}
all_filtered<-as.data.frame(do.call(rbind.data.frame, degs_filtered))
all_filtered_sig<-all_filtered[(all_filtered$padj<0.05),]
all_filtered_sig<-na.omit(all_filtered_sig)

all_filtered_sig<-all_filtered_sig%>% 
  mutate(Timepoint=case_when(
    comparison == "R1" ~ "24h",
    comparison == "R2" ~ "4d",
    comparison == "R3" ~ "8d",
    comparison == "R4" ~ "10d",
    comparison == "S1" ~ "24h",
    comparison == "S2" ~ "4d",
    comparison == "S3" ~ "8d",
    comparison == "S4" ~ "10d",
  ))

all_filtered_sig<-all_filtered_sig%>% 
  mutate(Genotype=case_when(
    comparison == "R1" ~ "Resistant",
    comparison == "R2" ~ "Resistant",
    comparison == "R3" ~ "Resistant",
    comparison == "R4" ~ "Resistant",
    comparison == "S1" ~ "Susceptible",
    comparison == "S2" ~ "Susceptible",
    comparison == "S3" ~ "Susceptible",
    comparison == "S4" ~ "Susceptible"
  ))

all_filtered_sig<-all_filtered_sig%>% 
  mutate(Regulation=case_when(
    log2FoldChange > 0 ~ "Upregulated",
    log2FoldChange < 0 ~ "Downregulated",
  ))

kable(table(all_filtered_sig$Genotype, 
            all_filtered_sig$Timepoint,
            all_filtered_sig$Regulation
            ))


degs<-as.character(unique(all_filtered_sig$row))
degs_all_data<-subset(all_filtered, all_filtered$row %in% degs)

degs_all_data_M<-degs_all_data[c(1, 3, 8)]
degs_all_data_M<-spread(degs_all_data_M, key="comparison", value="log2FoldChange", fill=0)

write.csv(all_filtered_sig, "~/Documents/Kat/Elites/data/separate_filtered.csv")
write.csv(degs_all_data_M, file="DEGs_all_data_matrix.csv")
write.csv(degs_all_data, file="DEGs_all_data.csv")

# ===========================
# do DE testing of everything together to compare
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor)
dds<-DESeq(dds)

RES_S1<-as.data.frame(results(dds, contrast=c("Factor","T1S","C1S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R1<-as.data.frame(results(dds, contrast=c("Factor","T1R","C1R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R2<-as.data.frame(results(dds, contrast=c("Factor","T2R","C2R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S2<-as.data.frame(results(dds, contrast=c("Factor","T2S","C2S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S3<-as.data.frame(results(dds, contrast=c("Factor","T3S","C3S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R3<-as.data.frame(results(dds, contrast=c("Factor","T3R","C3R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_S4<-as.data.frame(results(dds, contrast=c("Factor","T4S","C4S"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))
RES_R4<-as.data.frame(results(dds, contrast=c("Factor","T4R","C4R"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=TRUE))

RES_S1$Genotype<-"Susceptible"
RES_R1$Genotype<-"Resistant"
RES_G1$Genotype<-"Stigg"
RES_R2$Genotype<-"Resistant"
RES_S2$Genotype<-"Susceptible"
RES_S3$Genotype<-"Susceptible"
RES_R3$Genotype<-"Resistant"
RES_S4$Genotype<-"Susceptible"
RES_R4$Genotype<-"Resistant"

RES_R1$Timepoint<-1
RES_R2$Timepoint<-4
RES_R3$Timepoint<-8
RES_R4$Timepoint<-10
RES_S1$Timepoint<-1
RES_S2$Timepoint<-4
RES_S3$Timepoint<-8
RES_S4$Timepoint<-10

all<-rbind(RES_R1,
           RES_R2,
           RES_R3,
           RES_R4,
           RES_S1,
           RES_S2,
           RES_S3,
           RES_S4)
all_sig<-all[(all$padj<0.05),]
all_sig<-na.omit(all_sig)

all_sig<-all_sig%>% 
  mutate(Regulation=case_when(
    log2FoldChange > 0 ~ "Upregulated",
    log2FoldChange < 0 ~ "Downregulated",
  ))

kable(table(all_sig$Genotype, 
            all_sig$Timepoint,
            all_sig$Regulation
))

