options(width=200)

setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
load("MDS_PCA.Rdata")
load("DEGs.Rdata")
load("DEG.enrichments.Rdata")

tpms = read.table("10.09.2018.HNF1A_remap.gene.tpm.tsv",h=T,row.names=1)
tpm_ann = tpms[,1,drop=F]
tpms = tpms[,-1]


#############################################################################
#########      WGCNA analysis - each stage separately          ##############
#############################################################################

library("WGCNA")
allowWGCNAThreads()

####  DE_network

DE_tpms=tpms[,colnames(vstMat_DE)]
DE_use=rownames(DE_tpms)[which(rowSums(DE_tpms>=1)>=6)] # 17404 
DE_wgcna_vst=vstMat_DE[which(rownames(vstMat_DE) %in% DE_use),]

dim(DE_wgcna_vst)
# [1] 17404    18

### use sva ComBat function to regress out the batch effects of differentiations:
library(sva)
batch<-design_DE$diff
modcombat<-model.matrix(~1, data=design_DE)
DE_wgcna_vst= ComBat(dat=DE_wgcna_vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

degData = t(DE_wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

### pick power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "DE.WGCNA_softThreshold.pdf", width = 9, height = 5);
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

### power=3 would be fine...

net = blockwiseModules(degData, power = 3, maxBlockSize = 20000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       saveTOMs=T, verbose = 3)
table(net$colors)

# 0    1    2    3    4    5    6    7    8    9   10 
# 4077 5007 4819 1545  702  344  308  166  159  146  131 

save(net, file="net_DE.p3.M10.Rdata")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleLabels)$eigengenes
datME=moduleEigengenes(degData,moduleLabels)$eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

MEs = moduleEigengenes(degData, moduleLabels)$eigengenes

col= c(rep("darkred",12), rep("coral",6))
design_DE$HNF1A[order(design_DE$HNF1A)]
pdf("DE.module_eigengenes.barplots.pdf",width=10, height=3)
for (i in 1:length(unique(moduleColors))){
  which.module = unique(moduleLabels)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order(design_DE$HNF1A)], col=col, ylab="eigengene expression",xlab="array sample",main=which.module)
}
dev.off()

gene2module = data.frame(gene=colnames(degData), geneName=sapply(colnames(degData), function(x) as.character(gene_ann[x,])), module=moduleLabels)
write.table(gene2module, file = "DE.WGCNA.gene2module.txt", sep = "\t", row.names=F, quote=F)

sort(apply(datME, 2, function(x) t.test(x ~ design_DE$HNF1A)$p.value))
# ME1         ME2         ME7         ME3        ME10         ME0         ME9         ME8         ME5         ME6         ME4 
# 0.000218806 0.001158724 0.011270570 0.018687404 0.095045808 0.159150661 0.241819933 0.351340777 0.449570156 0.498222344 0.908125045 


gene2module$module=factor(gene2module$module)
lapply(HNF1A_DEGs, function(x) hypergeo_module(x, gene2module))
## module 2 contains most DEGs: p=1.10e-244

### PE network ### 
PE_tpms=tpms[,colnames(vstMat_PE)]

PE_use=rownames(PE_tpms)[which(rowSums(PE_tpms>=1)>=5)] # 19520 
PE_wgcna_vst=vstMat_PE[which(rownames(vstMat_PE) %in% PE_use),]

dim(PE_wgcna_vst)
# [1] 19520    17

### use sva ComBat function to regress out the batch effects of differentiations:
batch<-design_PE$diff
modcombat<-model.matrix(~1, data=design_PE)
PE_wgcna_vst= ComBat(dat=PE_wgcna_vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

degData = t(PE_wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

### pick power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "PE.WGCNA_softThreshold.pdf", width = 9, height = 5);
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

### power=6
net = blockwiseModules(degData, power = 6, maxBlockSize = 25000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       loadTOM =T, verbose = 3, deepSplit=1)
table(net$colors)

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30 
# 2246 2319 2115 1954 1068  977  943  818  759  717  695  679  464  434  369  364  338  312  269  254  219  192  155  151  149  128  104  100   87   77   64 
save(net, file="net_PE.p6.M30.Rdata")


table(moduleColors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleLabels)$eigengenes
datME=moduleEigengenes(degData,moduleLabels)$eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

MEs = moduleEigengenes(degData, moduleLabels)$eigengenes

col= c(rep("darkblue",12), rep("cornflowerblue",5))
design_PE$HNF1A[order(design_PE$HNF1A)]
pdf("PE.module_eigengenes.barplots.pdf",width=10, height=3)
for (i in 1:length(unique(moduleLabels))){
  which.module = unique(moduleLabels)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order(design_PE$HNF1A)], col=col, ylab="eigengene expression",xlab="array sample",main=which.module)
}
dev.off()

gene2module = data.frame(gene=colnames(degData), geneName=sapply(colnames(degData), function(x) as.character(gene_ann[x,])), module=moduleLabels)
write.table(gene2module, file = "PE.WGCNA.gene2module.txt", sep = "\t", row.names=F, quote=F)

sort(apply(datME, 2, function(x) t.test(x ~ design_PE$HNF1A)$p.value))
# ME5         ME19          ME3          ME9         ME12          ME2         ME29          ME8         ME30         ME20         ME26         ME14         ME10         ME23         ME11 
# 0.0001091669 0.0001314558 0.0001753764 0.0010803867 0.0017419022 0.0037992250 0.0098363071 0.0154593370 0.0178480342 0.0243008735 0.0342466515 0.0395868491 0.0696083678 0.0786378803 0.0865027303 
# ME22         ME25         ME24         ME15         ME17         ME13          ME4          ME0         ME21         ME16          ME7         ME28         ME27          ME1          ME6 
# 0.0988463956 0.1148433638 0.1475174870 0.1574000735 0.2001978376 0.2252731020 0.2368185927 0.2794696209 0.3387558438 0.3807489345 0.3827147456 0.4411730509 0.5477484928 0.6984423145 0.7427954838 
# ME18 
# 0.9943811265 
gene2module[gene2module$geneName=="HNF1A",] ## module:3
gene2module[gene2module$geneName=="PTF1A",] # M3
## M3 also contains PTF1A, NKX6-1, ONECUT1
gene2module[gene2module$geneName=="SOX9",] # M12

gene2module$module=factor(gene2module$module)
lapply(HNF1A_DEGs, function(x) hypergeo_module(x, gene2module))
## module 3 is enriched in PE up DEGs, also M12, M9
## PE down DEGs - enriched in M2, M5, M19
lapply(HNF1A_targets, function(x) hypergeo_module(x, gene2module))

######### EN network ################## 
EN_tpms=tpms[,colnames(vstMat_EN)]

EN_use=rownames(EN_tpms)[which(rowSums(EN_tpms>=1)>=6)] # 20327 
EN_wgcna_vst=vstMat_EN[which(rownames(vstMat_EN) %in% EN_use),]

dim(EN_wgcna_vst)
# [1] 20327    18

### use sva ComBat function to regress out the batch effects of differentiations:
batch<-design_EN$diff
modcombat<-model.matrix(~1, data=design_EN)
EN_wgcna_vst= ComBat(dat=EN_wgcna_vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

degData = t(EN_wgcna_vst)
gsg = goodSamplesGenes(degData, verbose = 3);
gsg$allOK

### pick power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
pdf(file = "EN.WGCNA_softThreshold.pdf", width = 9, height = 5);
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


net = blockwiseModules(degData, power = 5, maxBlockSize = 25000,
                       TOMType = "signed", networkType="signed",
                       numericLabels = TRUE, pamStage = T, pamRespectsDendro = T,
                       saveTOMs=T, verbose = 3)
table(net$colors)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 1722 3173 2149 2131 2070 1208  903  819  785  707  544  460  426  416  341  277  262  256  240  203  181  167  163  121  114   96   92   82   71   56   48   44 
save(net, file="net_EN.p5.M31.Rdata")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

lsdatME=moduleEigengenes(degData,moduleLabels)$eigengenes
datME=moduleEigengenes(degData,moduleLabels)$eigengenes
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )

MEs = moduleEigengenes(degData, moduleLabels)$eigengenes

col= c(rep("darkgreen",12), rep("green",6))
design_EN$HNF1A[order(design_EN$HNF1A)]
pdf("EN.module_eigengenes.barplots.pdf",width=10, height=3)
for (i in 1:length(unique(moduleLabels))){
  which.module = unique(moduleLabels)[i]
  ME=datME[, paste("ME",which.module, sep="")]
  barplot(ME[order(design_EN$HNF1A)], col=col, ylab="eigengene expression",xlab="array sample",main=which.module)
}
dev.off()

gene2module = data.frame(gene=colnames(degData), geneName=sapply(colnames(degData), function(x) as.character(gene_ann[x,])), module=moduleLabels)
write.table(gene2module, file = "EN.WGCNA.gene2module.txt", sep = "\t", row.names=F, quote=F)

sort(apply(datME, 2, function(x) t.test(x ~ design_EN$HNF1A)$p.value))
# ME5          ME1          ME3          ME4         ME18          ME2         ME11         ME12         ME17          ME9         ME25          ME8         ME28          ME6 
# 0.0004446998 0.0006558086 0.0027570902 0.0147703579 0.0332073088 0.0467153588 0.0710932189 0.1623016799 0.1803232936 0.1816875822 0.2460080045 0.2811047079 0.2852713570 0.3567586088 

gene2module[gene2module$geneName=="HNF1A",] ## module:9
gene2module[gene2module$geneName=="GCG",] # M1


## M1: LOXL4, RFX6, MAFB, GC, PTGER3, MTNR1B, SLC30A8, KCNJ11, CHGA, WFS1, NEUROD1, PAX4
## M2: FOXA2, NKX6-1, PDX1, 
## M4: ARX, PPY, GHRL
## M5: HNF1B, IAPP,PAX6, delta cell markers: SST, HHEX, CCKBR
## M9: HNF1A, MAFA, INS

gene2module$module=factor(gene2module$module)
lapply(HNF1A_DEGs, function(x) hypergeo_module(x, gene2module))
## EN up DEGs: M4
## EN down DEGs: M3, M5
lapply(HNF1A_targets, function(x) hypergeo_module(x, gene2module))
### M2 & M7  enriched  in HNF1A targets
### EndoC down (genes up-regulated by HNF1A): M1, M7, M9, M13
### EndoC up (genes down-regulated by HNF1A): M3, M6, M8, M10, M15, M17, M18, M19, 
### insulin secretion: M7, M9, M17
### alpha cells: M1
## beta: M9
### delta cells: M3, M5
## pp: M5, M8

save.image("HNF1A.WGCNA.Rdata")

options(stringsAsFactors = FALSE);
library("anRichment");

moduleLabels = gene2module$module
table(moduleLabels)
# Convert symbols to Entrez IDs
entrez = convert2entrez(organism = "human", symbol = gene2module$geneName);
# How many conversions were successful?
#table(is.finite(entrez))

GOcollection = buildGOcollection(organism = "human")
biosysCollection = BioSystemsCollection("human")

combinedCollection = mergeCollections(
  GOcollection,
  biosysCollection)

EN_GOenrichment = enrichmentAnalysis(
  classLabels = moduleLabels, identifiers = entrez,
  refCollection = combinedCollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = 0);

collectGarbage()

table.display = EN_GOenrichment$enrichmentTable;
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 70,
                                            split = "|");
head(table.display)
write.csv(EN_GOenrichment$enrichmentTable, file = "EN.GO_pathway_enrichment.anRichment.csv",
          row.names = FALSE)
