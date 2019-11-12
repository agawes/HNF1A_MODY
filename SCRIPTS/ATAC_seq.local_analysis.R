# grep -v ^# merged_peaks.counts.110918.txt | perl -i -pe 's/..\/merged_across_runs\/bams\///g' | perl -i -pe 's/.bam//g' > merged_peaks.optimal_set.counts.clean.110918.txt
# module load R/3.2.5

# R

## load libraries:
library(edgeR)
library("DESeq2")
library(sva)

setwd("~/Desktop/StemBANCC/HNF1A/ATAC-seq/")
counts = read.table("merged_peaks.optimal_set.counts.clean.110918.txt",h=T,row.names=1)
peak_ann = counts[,1:5,drop=F]
counts = counts[,-c(1:5)]

rename=read.table("design_table.txt",sep="\t",h=T)

## rename samples
colnames(counts) = as.character(sapply(colnames(counts), function(x) rename$Name[which(rename$Sample.ID == x)]))

### MDS plot
clone=sapply(strsplit(colnames(counts),split="_"), function(x) x[1])
diff=sapply(strsplit(colnames(counts),split="_"), function(x) x[2])
stage=sapply(strsplit(colnames(counts),split="_"), function(x) x[3])
hnf1a=rep("MODY",54)
hnf1a[grep("42-22|92-18",clone)]<-"corr"

design = data.frame(row.names = colnames(counts), clone=clone, diff = diff, stage=stage, HNF1A=hnf1a)



### MDS plot
design = design[colnames(counts),]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = design,
  design = ~1)

vsd <- varianceStabilizingTransformation(dds) 
vstMat = assay(vsd)

### col by stage and mutation - darker --> corrected
### DE - darkred & coral, PE - darkblue & cornflowerblue, EN - darkgreen & green
### pch by diff: 15,16,17

design$stage_mody = factor(paste(design$stage, design$HNF1A,sep="_"))
design$stage_mody = factor(design$stage_mody,levels(design$stage_mody)[c(1:2,5:6,3:4)])

col_stage = c("darkred","coral","darkblue","cornflowerblue","darkgreen","green")
col = sapply(design$stage_mody, function(x) col_stage[which(levels(design$stage_mody) == x)])

pch_diff=c(15,16,17)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])


### batch correction for pooling:
rename$batch=sapply(strsplit(as.character(rename$Sample.ID), split="_"), function(x) x[2])
design$pool=factor(sapply(rownames(design), function(x) rename[which(rename$Name==x),]$batch))

batch<-design$pool
modcombat<-model.matrix(~1, data=design)
combat_mydata= ComBat(dat=vstMat, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

pdf("MDS_all.pdf")
plotMDS(vstMat, col = col, pch = pch, cex=2)
legend(0.4,1.2,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Stage_MODY")
legend(1.0,1.2,levels(design$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat, col = col, pch = pch, cex=2, gene.selection="common")
legend(0.5,-0.6,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Protocol")
legend(1.25,-0.6,levels(design$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(combat_mydata, col = col, pch = pch, cex=2, main="Batch corrected")
legend(0.4,1.2,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Protocol")
legend(1.0,1.2,levels(design$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(combat_mydata, col = col, pch = pch, cex=2, gene.selection="common", main="Batch corrected")
legend(0.5,-0.6,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Protocol")
legend(1.25,-0.6,levels(design$diff),pch=pch_diff,bty="n",title="Diff")

dev.off()

### batch correction for differentiation:
batch<-design$diff
modcombat<-model.matrix(~1, data=design)
combat_mydata= ComBat(dat=combat_mydata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

pdf("PCA_all.sva_differentiation.pdf")
plotMDS(combat_mydata, col = col, pch = pch, cex=2, main="Diff corrected")
plotMDS(combat_mydata, col = col, pch = pch, cex=2, gene.selection="common", main="Diff corrected")
dev.off()


##### PCA's by stage:
col_clone = c("orange","blue","coral","red","cornflowerblue","darkred")
col = sapply(design$clone, function(x) col_clone[which(levels(design$clone) == x)])
col_order=c(6,1,3,4,2,5)

pch_diff=c(15,16,17)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])



design_DE=design[design$stage=="DE",]
dds_DE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="DE")],
  colData = design_DE,
  design = ~1)

vsd_DE <- varianceStabilizingTransformation(dds_DE) 
vstMat_DE = assay(vsd_DE)


pdf("MDS_DE.pdf")
par(mar=c(5.1,5.1,1,1))
# plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2)
# legend(0.4,0.45,levels(design_DE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
# legend(0.4,0.2,levels(design_DE$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2, gene.selection="common")
legend(0.4,0.4,levels(design_DE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.4,0.2,levels(design_DE$diff),pch=pch_diff,bty="n",title="Diff")

batch<-design_DE$diff
modcombat<-model.matrix(~1, data=design_DE)
cv=apply(vstMat_DE, 1, function(x) sd(x)/mean(x))
vstMat_DE_diff= ComBat(dat=vstMat_DE[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

x <- as.matrix(vstMat_DE_diff)
nsamples <- ncol(x)
cn <- colnames(x)
bad <- rowSums(is.finite(x)) < nsamples
nprobes <- nrow(x)
top <- min(500, nprobes)
labels <- colnames(x)
labels <- as.character(labels)
gene.selection <-  "common"
dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, 
                                                                  cn))
topindex <- nprobes - top + 1
for (i in 2:(nsamples)) for (j in 1:(i - 1)) dd[i, j] = sqrt(mean(sort.int((x[, 
                                                                              i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
a1 <- cmdscale(as.dist(dd), k = 2,eig=T)
total <- sum(abs(a1$eig))
100*sum(abs(a1$eig[1]))/total	    ## [1] 52.76693
100*sum(abs(a1$eig[2]))/total		## [1] 10.4133

plotMDS(vstMat_DE_diff, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2.5,
        gene.selection="common", cex.lab=2.5, xlab="PC1 (9.4 %)", ylab="PC2 (7.6 %)")
dev.off()



design_PE=design[design$stage=="PE",]
dds_PE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="PE")],
  colData = design_PE,
  design = ~1)

vsd_PE <- varianceStabilizingTransformation(dds_PE) 
vstMat_PE = assay(vsd_PE)


pdf("MDS_PE.pdf")
par(mar=c(5.1,5.1,1,1))

# plotMDS(vstMat_PE, col = col[which(design$stage=="PE")], pch = pch[which(design$stage=="PE")], cex=2)
# legend(-2.8,1.75,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
# legend(-2.8,0.45,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_PE, col = col[which(design$stage=="PE")], pch = pch[which(design$stage=="PE")], cex=2, gene.selection="common")
# legend(-2.8,1.75,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
# legend(-2.8,0.45,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")

batch<-design_PE$diff
modcombat<-model.matrix(~1, data=design_PE)
cv=apply(vstMat_PE, 1, function(x) sd(x)/mean(x))
vstMat_PE_diff= ComBat(dat=vstMat_PE[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

############### get variance explained by 1st and 2nd MDS dimension ############
x <- as.matrix(vstMat_PE_diff)
nsamples <- ncol(x)
cn <- colnames(x)
bad <- rowSums(is.finite(x)) < nsamples
nprobes <- nrow(x)
top <- min(500, nprobes)
labels <- colnames(x)
labels <- as.character(labels)
gene.selection <-  "common"
dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, 
                                                                  cn))
topindex <- nprobes - top + 1
for (i in 2:(nsamples)) for (j in 1:(i - 1)) dd[i, j] = sqrt(mean(sort.int((x[, 
                                                                              i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
a1 <- cmdscale(as.dist(dd), k = 2,eig=T)
total <- sum(abs(a1$eig))
100*sum(abs(a1$eig[1]))/total	    ## [1] 52.76693
100*sum(abs(a1$eig[2]))/total		## [1] 10.4133


plotMDS(vstMat_PE_diff, col = col[which(design$stage=="PE")], pch = pch[which(design$stage=="PE")], cex=2.5, 
        gene.selection="common", cex.lab=2.5, xlab="PC1 (9.8%)", ylab="PC2 (8.3%)")
dev.off()


design_EN=design[design$stage=="EN",]
dds_EN <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="EN")],
  colData = design_EN,
  design = ~1)

vsd_EN <- varianceStabilizingTransformation(dds_EN) 
vstMat_EN = assay(vsd_EN)

pdf("MDS_EN.pdf")
par(mar=c(5.1,5.1,1,1))

# plotMDS(vstMat_EN, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2)
# legend(0.25,0.7,levels(design_EN$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
# legend(0.25,0.35,levels(design_EN$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_EN, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2, gene.selection="common")
legend(0.25,-0.1,levels(design_EN$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.25,-0.45,levels(design_EN$diff),pch=pch_diff,bty="n",title="Diff")

batch<-design_EN$diff
modcombat<-model.matrix(~1, data=design_EN)
cv=apply(vstMat_EN, 1, function(x) sd(x)/mean(x))
vstMat_EN_diff= ComBat(dat=vstMat_EN[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

x <- as.matrix(vstMat_EN_diff)
nsamples <- ncol(x)
cn <- colnames(x)
bad <- rowSums(is.finite(x)) < nsamples
nprobes <- nrow(x)
top <- min(500, nprobes)
labels <- colnames(x)
labels <- as.character(labels)
gene.selection <-  "common"
dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, 
                                                                  cn))
topindex <- nprobes - top + 1
for (i in 2:(nsamples)) for (j in 1:(i - 1)) dd[i, j] = sqrt(mean(sort.int((x[, 
                                                                              i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
a1 <- cmdscale(as.dist(dd), k = 2,eig=T)
total <- sum(abs(a1$eig))
100*sum(abs(a1$eig[1]))/total	    ## [1] 52.76693
100*sum(abs(a1$eig[2]))/total		## [1] 10.4133

plotMDS(vstMat_EN_diff, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2.5, 
        gene.selection="common", cex.lab=2.5, xlab="PC1 (9.6 %)", ylab="PC2 (7.5 %)")

dev.off()

save.image("hnf1a.atac.Rdata")


############ differential openness between MODY and control  #############
############ 	separately for each stage			 #############
### adjust for diff and pool ###

################ DE ################

design_DE = design[grep("DE",design$stage),]
design_DE$HNF1A=factor(design_DE$HNF1A, levels=c("MODY","corr"))

dds_DE <- DESeqDataSetFromMatrix(
  countData = counts[,grep("DE",design$stage)],
  colData = design_DE,
  design = ~pool+diff+HNF1A)

dds_DE <- DESeq(dds_DE)
res_DE = results(dds_DE)
length(which(res_DE$padj < 0.05))	### 27
length(which(res_DE$padj < 0.05 & res_DE$log2FoldChange>0))	### 9
length(which(res_DE$padj < 0.05 & res_DE$log2FoldChange<0))	### 18

res_DE = cbind(peak_ann, res_DE)
head(data.frame(res_DE[order(res_DE$padj),]),n=20)
write.table(res_DE, file = "DE.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

# top peak: peak_chr12_126467410; padj=2.052324e-35 --> next p=1.748439e-09

################ PE ################

design_PE = design[grep("PE",design$stage),]
design_PE$HNF1A=factor(design_PE$HNF1A, levels=c("MODY","corr"))

dds_PE <- DESeqDataSetFromMatrix(
  countData = counts[,grep("PE",design$stage)],
  colData = design_PE,
  design = ~pool+diff+HNF1A)

dds_PE <- DESeq(dds_PE)
res_PE = results(dds_PE)
length(which(res_PE$padj < 0.05))	### 7811
length(which(res_PE$padj < 0.05 & res_PE$log2FoldChange>0))	### 7157
length(which(res_PE$padj < 0.05 & res_PE$log2FoldChange<0))	### 654

res_PE = cbind(peak_ann, res_PE)
head(data.frame(res_PE[order(res_PE$padj),]),n=20)
write.table(res_PE, file = "PE.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## second best peak:  peak_chr12_126467410
## HNF1A: chr12:121,416,346-121,440,315

################ EN ################

design_EN = design[grep("EN",design$stage),]
design_EN$HNF1A=factor(design_EN$HNF1A, levels=c("MODY","corr"))

dds_EN <- DESeqDataSetFromMatrix(
  countData = counts[,grep("EN",design$stage)],
  colData = design_EN,
  design = ~pool+diff+HNF1A)

dds_EN <- DESeq(dds_EN)
res_EN = results(dds_EN)
length(which(res_EN$padj < 0.05))	### 1710
length(which(res_EN$padj < 0.05 & res_EN$log2FoldChange>0))	### 1681
length(which(res_EN$padj < 0.05 & res_EN$log2FoldChange<0))	### 29

res_EN = cbind(peak_ann, res_EN)
head(data.frame(res_EN[order(res_EN$padj),]),n=20)
write.table(res_EN, file = "EN.DOCs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

## top peak: peak_chr12_126467410 
peak_chr12_126467410 chr12 126467410 126469059      +   1650

###############################################################

# boxplot of the top peak depth:

pdf("top_peak.box.pdf", height=4)
boxplot(log10(as.numeric(counts["peak_chr12_126467410",])) ~ design$stage_mody, col=col_stage, ylab = "log10(counts)", main="peak_chr12_126467410")
dev.off()





#######################################


#### run WGCNA
library("WGCNA")
allowWGCNAThreads()

### only peaks with at least 30x in >= 6 samples & with CV > 0.5
peak_CV = apply(counts, 1, function(x) sd(x)/mean(x))
length(intersect(which(peak_CV>0.5),which(rowSums(counts>=30)>=6)))
[1] 30542

use=rownames(counts)[intersect(which(peak_CV>0.5),which(rowSums(counts>=30)>=6))] # 105621
wgcna_vst=vstMat[which(rownames(vstMat) %in% use),]
dim(wgcna_vst)
# [1] 30542    54

degData = t(wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

### pick power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "WGCNA_softThreshold.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()




##### run WGCNA module detection

net = blockwiseModules(degData, power = 12, maxBlockSize = 40000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       saveTOMs=T, verbose = 3)
table(net$colors)

0     1     2     3     4     5
36 13668  9850  3833  2660   495




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleColors)$eigengenes
datME=moduleEigengenes(degData,moduleColors)$eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

MEs = moduleEigengenes(degData, moduleColors)$eigengenes
order=order(design$stage_mody, design$clone, design$diff)

pdf("Module_eigengenes.barplots.pdf",width=10, height=4)
par(las=2, mar=c(8,5,4,1))
for (i in 1:length(unique(moduleColors))){
  which.module = unique(moduleColors)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order], col=col[order], ylab="eigengene expression",xlab="",main=which.module, names.arg=rownames(design)[order], cex.names=0.75)
}
dev.off()

### as in RNA-seq the detected modules represent stage-specific variation
### to dissect the MODY-WT differences, we need to run 3 separate networks for each stage



peak2module = data.frame(peak=colnames(degData), peak_genes=sapply(colnames(degData), function(x) as.character(peak_ann[x,])), module=moduleColors)
write.table(peak2module, file = "WGCNA.peak2module.all_stages.txt", sep = "\t", row.names=F, quote=F)





######### correlate the RNA and ATAC modules #########

rownames(datME) = colnames(wgcna_vst)
write.table(datME, file="ATAC.module_eigengenes.endo_comp.txt",sep="\t",quote=F)
rna_datME=read.table("../SmartSeq/RNA.module_eigengenes.endo_comp.txt",h=T,row.names=1)
rownames(rna_datME) = gsub("GFP\\.","GFP-", rownames(rna_datME))
rna_datME= rna_datME[rownames(datME),]
write.table(cor(datME, rna_datME), file="RNA_vs_ATAC.module_eigengene_cor.txt",sep="\t",quote=F)


cormat=round(cor(datME, rna_datME),2)

#cluster cor's
hc=heatmap(cormat)
cormat=cormat[hc$rowInd,hc$colInd]

library(reshape2)
library(ggplot2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

melted_cormat <- melt(cormat)

# Heatmap
pdf("ATAC_RNA_cor_matrix.pdf")
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +xlab("RNA-seq") + ylab("ATAC-seq")

dev.off()

save.image("atac.wgcna.Rdata")

