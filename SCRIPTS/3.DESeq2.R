## load libraries:
library(edgeR)
library("DESeq2")
options(width=200)


setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
counts = read.table("10.09.2018.HNF1A_remap.gene.counts.tsv",h=T,row.names=1)
gene_ann = counts[,1,drop=F]
counts = counts[,-1]

clone=sapply(strsplit(colnames(counts),split="_"), function(x) x[1])
diff=sapply(strsplit(colnames(counts),split="_"), function(x) x[2])
stage=sapply(strsplit(colnames(counts),split="_"), function(x) x[3])
hnf1a=c(rep("corrected",9), rep("MODY",27), rep("corrected",9), rep("MODY",9))

design = data.frame(row.names = colnames(counts), clone=clone, diff = diff, stage=stage, HNF1A=hnf1a)


#############################################################################
#################    Differential expression: DESeq2        #################
#############################################################################

design_DE=design[design$stage=="DE",]
design_DE$HNF1A=factor(design_DE$HNF1A, levels=c("MODY","corrected"))
dds_DE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="DE")],
  colData = design_DE,
  design = ~diff+HNF1A)

dds_DE <- DESeq(dds_DE)
res_DE = results(dds_DE)
length(which(res_DE$padj < 0.05))	### 1841
length(which(res_DE$padj < 0.05 & res_DE$log2FoldChange>0)) ## 610
length(which(res_DE$padj < 0.01 & abs(res_DE$log2FoldChange)>=1)) ## 57

res_DE = cbind(gene_ann, res_DE)

head(data.frame(res_DE[order(res_DE$padj),]),n=20)
write.table(res_DE, file = "HNF1A_DE.DEGs.DESeq2.txt",sep="\t",quote=F,col.names=NA)


#########  PE #########

design_PE=design[design$stage=="PE",][-3,]
design_PE$HNF1A=factor(design_PE$HNF1A, levels=c("MODY","corrected"))

dds_PE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="PE")][,-3],
  colData = design_PE,
  design = ~diff+HNF1A)

dds_PE <- DESeq(dds_PE)
res_PE = results(dds_PE)
length(which(res_PE$padj < 0.05))	### 340
length(which(res_PE$padj < 0.05 & res_PE$log2FoldChange>0)) ## 169
length(which(res_PE$padj < 0.01 & abs(res_PE$log2FoldChange)>=1)) ## 11

res_PE = cbind(gene_ann, res_PE)

head(data.frame(res_PE[order(res_PE$padj),]),n=20)
write.table(res_PE, file = "HNF1A_PE.DEGs.DESeq2.txt",sep="\t",quote=F,col.names=NA)


#### EN #####
design_EN=design[design$stage=="EN",]
design_EN$HNF1A=factor(design_EN$HNF1A, levels=c("MODY","corrected"))

dds_EN <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="EN")],
  colData = design_EN,
  design = ~diff+HNF1A)

dds_EN <- DESeq(dds_EN)
res_EN = results(dds_EN)
length(which(res_EN$padj < 0.05))	### 1870
length(which(res_EN$padj < 0.05 & res_EN$log2FoldChange>0)) ## 949
length(which(res_EN$padj < 0.01 & abs(res_EN$log2FoldChange)>=1)) ## 132

res_EN = cbind(gene_ann, res_EN)

head(data.frame(res_EN[order(res_EN$padj),]),n=20)
write.table(res_EN, file = "HNF1A_EN.DEGs.DESeq2.txt",sep="\t",quote=F,col.names=NA)

save.image("DEGs.Rdata")