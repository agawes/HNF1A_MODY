library(beeswarm)

setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
counts = read.table("10.09.2018.HNF1A_remap.gene.counts.tsv",h=T,row.names=1)
gene_ann = counts[,1,drop=F]
counts = counts[,-1]

tpms = read.table("10.09.2018.HNF1A_remap.gene.tpm.tsv",h=T,row.names=1)
tpm_ann = tpms[,1,drop=F]
tpms = tpms[,-1]

clone=sapply(strsplit(colnames(counts),split="_"), function(x) x[1])
diff=sapply(strsplit(colnames(counts),split="_"), function(x) x[2])
stage=sapply(strsplit(colnames(counts),split="_"), function(x) x[3])
hnf1a=c(rep("corrected",9), rep("MODY",27), rep("corrected",9), rep("MODY",9))

design = data.frame(row.names = colnames(counts), clone=clone, diff = diff, stage=stage, HNF1A=hnf1a)

design$stage_mody = factor(paste(design$stage, gsub("ected","",design$HNF1A),sep="_"))
design$stage_mody = factor(design$stage_mody,levels(design$stage_mody)[c(1:2,5:6,3:4)])

col_stage = c("darkred","coral","darkblue","cornflowerblue","darkgreen","green")
col = sapply(design$stage_mody, function(x) col_stage[which(levels(design$stage_mody) == x)])

pch_diff=c(15,16,17)
pch = sapply(design$diff, function(x) pch_diff[which(levels(design$diff) == x)])

#############################################################################
#####################      Gene expression plots         ####################
#############################################################################

## reorder & rename the factor levels for plotting
design$stage_mody = factor(design$stage_mody, levels=c("DE_MODY","DE_corr","PE_MODY","PE_corr","EN_MODY","EN_corr"))
levels(design$stage_mody)=c("DE_MODY","DE_WT_GE","PE_MODY","PE_WT_GE","EN_MODY","EN_WT_GE")

gene="RBP4"
#pdf(paste0(gene,".pdf"),width=12)
par(mar=c(5.1,5.1,4,1))
### this will make a plot with different symbols for each rep:
beeswarm(as.numeric(tpms[which(tpm_ann$GeneName==gene),]) ~ design$stage_mody, main=gene, pwcol=col, pwpch=pch, 
         ylab = paste(gene,"TPMs"), xlab="", cex=2.5, cex.lab=2, cex.main=2, cex.axis=1.75, bty="l",xaxt="n")
axis(1,at=1:6,labels=sub("_","\n",levels(design$stage_mody)), cex.axis=2, padj=0.8)
bxplot(as.numeric(tpms[which(tpm_ann$GeneName==gene),]) ~ design$stage_mody, add = TRUE)
#dev.off()

### single stage plots
col=col[which(rownames(design) != "Clone42.22_Diff.012_PE")]
design1=design[which(rownames(design) != "Clone42.22_Diff.012_PE"),] ## rm the outlier
tpms=tpms[,which(rownames(design) != "Clone42.22_Diff.012_PE")]

stage="PE"
gene="CACNA1A"
tmp=as.numeric(tpms[which(tpm_ann$GeneName==gene),grepl(stage,design1$stage)])
group=factor(as.character(design1$stage_mody[grepl(stage,design1$stage)]))
pdf(paste0(gene,".pdf"),width=5,heigh=4)
par(mar=c(5.1,5.1,4,1))
### this will make a plot with different symbols for each rep:
beeswarm(tmp ~ group, main="", pwcol=col[grepl(stage,design1$stage)], pwpch=pch[grepl(stage,design1$stage)], 
         ylab = paste(gene,"TPMs"), xlab="", cex=2.5, cex.lab=2, cex.main=2, cex.axis=1.75, bty="l",xaxt="n")
axis(1,at=1:2,labels=sub("_","\n",levels(group)), cex.axis=2, padj=0.8)
bxplot(tmp ~ group, add = TRUE)
dev.off()

