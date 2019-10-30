## load libraries:
library(edgeR)
library("DESeq2")
options(width=200)

setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
counts = read.table("10.09.2018.HNF1A_remap.gene.counts.tsv",h=T,row.names=1)
gene_ann = counts[,1,drop=F]
counts = counts[,-1]

#############################################################################
#####################     Global PCA/MDS analysis        ####################
#############################################################################

clone=sapply(strsplit(colnames(counts),split="_"), function(x) x[1])
diff=sapply(strsplit(colnames(counts),split="_"), function(x) x[2])
stage=sapply(strsplit(colnames(counts),split="_"), function(x) x[3])
hnf1a=c(rep("corrected",9), rep("MODY",27), rep("corrected",9), rep("MODY",9))

design = data.frame(row.names = colnames(counts), clone=clone, diff = diff, stage=stage, HNF1A=hnf1a)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = design,
  design = ~1)

vsd <- varianceStabilizingTransformation(dds) 
vstMat = assay(vsd)

### col by stage and mutation - darker --> corrected
### DE - darkred & coral, PE - darkblue & cornflowerblue, EN - darkgreen & green
### pch by diff: 15,16,17

design$stage_mody = factor(paste(design$stage, gsub("ected","",design$HNF1A),sep="_"))
design$stage_mody = factor(design$stage_mody,levels(design$stage_mody)[c(1:2,5:6,3:4)])

col_stage = c("darkred","coral","darkblue","cornflowerblue","darkgreen","green")
col = sapply(design$stage_mody, function(x) col_stage[which(levels(design$stage_mody) == x)])

pch_diff=c(15,16,17)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])

pdf("MDS_all.pdf")
plotMDS(vstMat, col = col, pch = pch, cex=2)
legend(-2.8,2.25,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Stage/MODY")
legend(-2.8,0.45,levels(design$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat, col = col, pch = pch, cex=2, gene.selection="common")
legend(-2.8,1.75,levels(design$stage_mody),col=col_stage, pch=16,bty="n", title="Stage/MODY")
legend(-2.8,0.45,levels(design$diff),pch=pch_diff,bty="n",title="Diff")
dev.off()

### clusters by stage

### per stage MDS

### col by clone: reds --> MODY; blue --> corrected
### SCF012 - darkred; Clone4 - orange; Clone5 - coral, Clone6 - red
### Clone42.22 - blue, Clone92.18 - cornflowerblue
### pch by diff: 15,16,17

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
plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2)
legend(0.4,0.45,levels(design_DE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.4,0.2,levels(design_DE$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2, gene.selection="common")
legend(0.4,0.4,levels(design_DE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.4,0.2,levels(design_DE$diff),pch=pch_diff,bty="n",title="Diff")
dev.off()



design_PE=design[design$stage=="PE",]
dds_PE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="PE")],
  colData = design_PE,
  design = ~1)

vsd_PE <- varianceStabilizingTransformation(dds_PE) 
vstMat_PE = assay(vsd_PE)


pdf("MDS_PE.with_outlier.pdf")
plotMDS(vstMat_PE, col = col[which(design$stage=="PE")], pch = pch[which(design$stage=="PE")], cex=2)
legend(-2.8,1.75,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(-2.8,0.45,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_PE, col = col[which(design$stage=="PE")], pch = pch[which(design$stage=="PE")], cex=2, gene.selection="common")
legend(-2.8,1.75,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(-2.8,0.45,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")
dev.off()


#### one of the PE samples is an outlier: Clone42.22 PE - diff012
design_PE=design_PE[-3,]

dds_PE <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="PE")][,-3],
  colData = design_PE,
  design = ~1)

vsd_PE <- varianceStabilizingTransformation(dds_PE) 
vstMat_PE = assay(vsd_PE)

pdf("MDS_PE.pdf")
plotMDS(vstMat_PE, col = col[which(design$stage=="PE")][-3], pch = pch[which(design$stage=="PE")][-3], cex=2)
legend(-0.4,0.7,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(-0.4,0.27,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_PE, col = col[which(design$stage=="PE")][-3], pch = pch[which(design$stage=="PE")][-3], cex=2, gene.selection="common")
legend(-0.4,0.65,levels(design_PE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(-0.4,0.25,levels(design_PE$diff),pch=pch_diff,bty="n",title="Diff")
dev.off()



design_EN=design[design$stage=="EN",]
dds_EN <- DESeqDataSetFromMatrix(
  countData = counts[,which(design$stage=="EN")],
  colData = design_EN,
  design = ~1)

vsd_EN <- varianceStabilizingTransformation(dds_EN) 
vstMat_EN = assay(vsd_EN)

pdf("MDS_EN.pdf")
plotMDS(vstMat_EN, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2)
legend(0.25,0.7,levels(design_EN$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.25,0.35,levels(design_EN$diff),pch=pch_diff,bty="n",title="Diff")

plotMDS(vstMat_EN, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2, gene.selection="common")
legend(0.25,-0.1,levels(design_EN$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.25,-0.45,levels(design_EN$diff),pch=pch_diff,bty="n",title="Diff")
dev.off()

######   regress out the effects of separate differentiations #########

library(sva)
batch<-design_DE$diff
modcombat<-model.matrix(~1, data=design_DE)
cv=apply(vstMat_DE, 1, function(x) sd(x)/mean(x))
vstMat_DE= ComBat(dat=vstMat_DE[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

############### get variance explained by 1st and 2nd MDS dimension ############
x <- as.matrix(vstMat_DE)
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

pdf("PCA_DE.sva_differentiation.pdf")
par(mar=c(5.1,5.1,1,1))
plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2.5, gene.selection="common",
        xlab= "PC1 (52.7%)", ylab = "PC2 (10.4%)", cex.lab=2.5)
plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2, gene.selection="common")
legend(-0.2,-0.1,levels(design_DE$clone)[col_order],col=col_clone[col_order], pch=16,bty="n", title="Clone")
legend(0.2,-0.1,levels(design_DE$diff),pch=pch_diff,bty="n",title="Diff")
plotMDS(vstMat_DE, col = col[which(design$stage=="DE")], pch = pch[which(design$stage=="DE")], cex=2, gene.selection="common")
legend(-0.2,-0.1,c("Pro291fsinsC/+ patient line","Pro291fsinsC/+ clone1","Pro291fsinsC/+ clone2", "Pro291fsinsC/+ clone3", "Corrected/+ clone1",
                   "Corrected/+ clone2"), col=col_clone[col_order], pch=16,bty="n", title="Clone")

dev.off()

batch<-design_PE$diff
modcombat<-model.matrix(~1, data=design_PE)
cv=apply(vstMat_PE, 1, function(x) sd(x)/mean(x))
vstMat_PE= ComBat(dat=vstMat_PE[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

############### get variance explained by 1st and 2nd MDS dimension ############
x <- as.matrix(vstMat_PE)
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
pdf("PCA_PE.sva_differentiation.pdf")
par(mar=c(5.1,5.1,1,1))
plotMDS(vstMat_PE, col = col[which(design$stage=="PE")][-3], pch = pch[which(design$stage=="PE")][-3], cex=2.5, 
        gene.selection="common", xlab="PC1 (44.4%)", ylab="PC2 (16.9%)", cex.lab=2.5)
dev.off()

batch<-design_EN$diff
modcombat<-model.matrix(~1, data=design_EN)
cv=apply(vstMat_EN, 1, function(x) sd(x)/mean(x))
vstMat_EN= ComBat(dat=vstMat_EN[cv>0,], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

x <- as.matrix(vstMat_EN)
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

pdf("PCA_EN.sva_differentiation.pdf")
par(mar=c(5.1,5.1,1,1))
plotMDS(vstMat_EN, col = col[which(design$stage=="EN")], pch = pch[which(design$stage=="EN")], cex=2.5, gene.selection="common",
        xlab ="PC1 (43.0%)", ylab="PC2 (12.8%)", cex.lab=2.5)
dev.off()

save.image("MDS_PCA.Rdata")
