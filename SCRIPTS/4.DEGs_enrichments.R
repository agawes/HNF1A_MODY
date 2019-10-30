setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
load("DEGs.Rdata")

#############################################################################
######    Enrichment of HNF1A targets and cell type markers in DEGs    ######
#############################################################################

#### check if DEGs are enriched in HNF1A targets
#### - from EndoC knockdowns
#### - from public databases - compiled in the tftargets package: TRRUST and Marbach2016
#### TRRUST - manuallly curated
#### Marbach - gene regulatory networks, tissue specific - but in this package it seems to be compiled across all tissues/cell types
#### there are 2 relevant files: pancreas_adult.txt and embryonic_pancreas_cell_line.tx among the 394 networks downloaded

HNF1A_targets = list(
  TRRUST=TRRUST$HNF1A, 		# 26
  Marbach=Marbach2016$HNF1A,   # 2069
  Marbach_embryonic_pancreas=scan("external_data/HNF1A_targets.embryonic_pancreas_cell_line.Marbach2016.txt", what="character"), # 832
  Marbach_pancreas_adult=scan("external_data/HNF1A_targets.pancreas_adult.Marbach2016.txt", what="character"), # 719
  EndoC_up=scan("external_data/endoC_KD.HNF1A.genes_up.fdr0.05.txt", what="character"), # 790 
  EndoC_down=scan("external_data/endoC_KD.HNF1A.genes_down.fdr0.05.txt",what="character"), # 653
  ins_secretion=scan("external_data/GO_0030073.insulin_secretion.txt",what="character"), # 202
  alpha=scan("external_data/Muraro2016.alpha_cells.txt",what="character"),
  beta=scan("external_data/Muraro2016.beta_cell.txt",what="character"),
  delta=scan("external_data/Muraro2016.delta_cells.txt",what="character"),
  pp=scan("external_data/Muraro2016.PP_cells.txt",what="character"),
  PP_Cebola=scan("~/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18/external_data/TEAD_YAP.500genes.txt",what="character"),
  PP_our=scan("~/Desktop/StemBANCC/endodermal competence/FINAL_Reanalysis_October18/PE_iPSC_signature.GFPpos.txt",what="character")
)
## up - means up in corrected clones
HNF1A_DEGs = list(
  DE_up = as.character(res_DE$GeneName[which(res_DE$padj<0.05 & res_DE$log2FoldChange>0)]),
  DE_down = as.character(res_DE$GeneName[which(res_DE$padj<0.05 & res_DE$log2FoldChange<0)]),
  PE_up = as.character(res_PE$GeneName[which(res_PE$padj<0.05 & res_PE$log2FoldChange>0)]),
  PE_down = as.character(res_PE$GeneName[which(res_PE$padj<0.05 & res_PE$log2FoldChange<0)]),
  EN_up = as.character(res_EN$GeneName[which(res_EN$padj<0.05 & res_EN$log2FoldChange>0)]),
  EN_down = as.character(res_EN$GeneName[which(res_EN$padj<0.05 & res_EN$log2FoldChange<0)])
)

##############  enrichment of genes in a gene list within modules ################
## run hypergeometric enrichment of a gene vector within a set of WGCNA modules

hypergeo<-function(genes, deg_list){
  enrichment = rep(1, length(deg_list))
  for(i in 1:length(deg_list)) {
    A=length(intersect(genes, deg_list[[i]] )) ## overlap
    B=length(deg_list[[i]])## module size
    C=nrow(res_DE)
    D=length(genes)
    enrichment[i]= phyper(A-1,B,C-B,D,lower.tail=F)
  }
  names(enrichment) = names(deg_list)
  return(enrichment)
}

## run hypergeometric enrichment of a gene vector within a set of WGCNA modules
hypergeo_module<-function(genes, gene2module){
  enrichment = rep(1, length(levels(gene2module$module)))
  for(i in 1:length(unique(gene2module$module))) {
    A=length(intersect(genes, gene2module$geneName[gene2module$module==levels(gene2module$module)[i]] )) ## overlap
    B=length(gene2module$geneName[gene2module$module==levels(gene2module$module)[i]])## module size
    C=length(gene2module$geneName)
    D=length(which(genes %in% gene2module$geneName))
    enrichment[i]= phyper(A-1,B,C-B,D,lower.tail=F)
  }
  names(enrichment) = levels(gene2module$module)
  return(enrichment)
}

lapply(HNF1A_targets, function(x) hypergeo(x, HNF1A_DEGs))

pdf("HNF1A_EndoC_target_enrichment_in_DEGs.pdf", height=5)
par(mfrow=c(2,1), mar=c(2.5,5,0.5,0))
barplot(-log10(hypergeo(HNF1A_targets$EndoC_down, HNF1A_DEGs)), ylab="-log10 p-value", ylim=c(0,135))
barplot(-log10(hypergeo(HNF1A_targets$EndoC_up, HNF1A_DEGs)), ylab="-log10 p-value", ylim=c(0,135))
dev.off()

#### correlation of EndoC and MODY/Corrected fold-changes
endo_res=read.table("external_data/EndoC.HNF1A_vs_NT.DEGs.DESeq2.tsv",sep="\t",h=T,row.names=1)
res_EN=data.frame(res_EN)
merged=merge(res_EN, endo_res,by="row.names",suffixes=c(".iPSC",".EndoC"))

pdf("endoC_target_overlap.pdf")
col=rep(rgb(0,0,0,50,maxColorValue = 255),nrow(merged))
col[which(merged$padj.EndoC<0.05 & merged$padj.iPSC<0.05)]="red"
par(mar=c(5.1,8.1,1,1))
plot(merged$log2FoldChange.EndoC,merged$log2FoldChange.iPSC,pch=16, col=col, xlab="log2FC: EndoC KD vs NT control",
     ylab="log2FC: Corr/+ vs Pro291fsInsC/+\n(iPSC model: EN st7)", cex.lab=2, cex.axis=2, ylim=c(-5,3))
dev.off()
#############################################################################
#########    Twister plots - insulin secretion and others        ##############
###########################################################################

#### twister plot of insulin secretion pathway:
intersect(HNF1A_targets$ins_secretion, HNF1A_DEGs$EN_up)
# [1] "AACS"    "BMP8A"   "CACNA1A" "CASR"    "CHGA"    "DPP4"    "GAL"     "GCG"     "GPR119"  "GPR27"   "HMGCR"
# [12] "HNF1A"   "ITPR1"   "KCNB1"   "KCNJ11"  "MPC2"    "MTNR1B"  "OXCT1"   "PTPRN"   "PTPRN2"  "RAPGEF4" "RFX6"
# [23] "SLC30A8" "SREBF1"
intersect(HNF1A_targets$ins_secretion, HNF1A_DEGs$EN_down)
# [1] "ADCYAP1" "CRH"     "ENSA"    "FAM3B"   "HNF1B"   "INHBB"   "NOS2"    "PIM3"    "REST"    "SYTL4"

plot_ins_secr_genes = c(intersect(HNF1A_targets$ins_secretion, HNF1A_DEGs$EN_up),intersect(HNF1A_targets$ins_secretion, HNF1A_DEGs$EN_down))

a=data.frame(res_EN[which(res_EN$GeneName %in% plot_ins_secr_genes),c(1,3,7)])
a=a[order(a$log2FoldChange),]
a$col=rep("grey",nrow(a))
a$col[which(a$log2FoldChange>0)]="darkred"
a$col[which(a$log2FoldChange<0)]="darkblue"

# a$hjust <- ifelse(a$log2FoldChange > 0, -0.5, 0.5)
# pdf("twister.HNF1A_EN.insulin_secretion.pdf",width=16)
# par(las=2)
# bar=barplot(a$log2FoldChange, ylab="log2FoldChange", col=a$col)
# text(bar,a$hjust, a$GeneName, srt=90, cex=1.3)
# dev.off()

a$hjust <- ifelse(a$log2FoldChange > 0, -0.5, 0.5)
pdf("twister.HNF1A_EN.insulin_secretion.pdf",height=16)
bar=barplot(a$log2FoldChange, xlab="log2FoldChange", col=a$col, horiz=T)
text(a$hjust,bar, a$GeneName, cex=1.3)
dev.off()

### twister plot of alpha cells genes:
a=data.frame(res_EN[which(res_EN$GeneName %in% HNF1A_targets$alpha),c(1,3,7)])
a=a[order(a$log2FoldChange),]
a$col=rep("grey",nrow(a))
a$col[which(a$padj<0.05 & a$log2FoldChange>0)]="darkred"
a$col[which(a$padj<0.05 & a$log2FoldChange<0)]="darkblue"

a$hjust <- ifelse(a$log2FoldChange > 0, -0.25, 0.25)
pdf("twister.HNF1A_EN.alpha.pdf",height=10)
par(mar=c(5.1,8,1,1), xpd=T)
bar=barplot(a$log2FoldChange, xlab="log2FoldChange", col=a$col, cex.lab=1.5, cex.axis=1.5, horiz=T)
text(a$hjust,bar, a$GeneName, cex=1.5)
dev.off()

plot_PP_genes = c(intersect(HNF1A_targets$PP_Cebola, HNF1A_DEGs$PE_up),intersect(HNF1A_targets$PP_Cebola, HNF1A_DEGs$PE_down))
a=data.frame(res_PE[which(res_PE$GeneName %in% plot_PP_genes),c(1,3,7)])
a=a[order(a$log2FoldChange),]
a$col=rep("grey",nrow(a))
a$col[which(a$padj<0.05 & a$log2FoldChange>0)]="darkred"
a$col[which(a$padj<0.05 & a$log2FoldChange<0)]="darkblue"

a$hjust <- ifelse(a$log2FoldChange > 0, -0.15, 0.15)
pdf("twister.HNF1A_PE.Cebola_PP_signature.pdf",height=10)
par(mar=c(5.1,8,1,1), xpd=T)
bar=barplot(a$log2FoldChange, xlab="log2FoldChange", col=a$col, cex.lab=1.5, cex.axis=1.5, horiz=T)
text(a$hjust,bar, a$GeneName, cex=1.5)
dev.off()

plot_PP_genes_2 = c(intersect(HNF1A_targets$PP_our, HNF1A_DEGs$PE_up),intersect(HNF1A_targets$PP_our, HNF1A_DEGs$PE_down))
a=data.frame(res_PE[which(res_PE$GeneName %in% plot_PP_genes_2),c(1,3,7)])
a=a[order(a$log2FoldChange),]
a$col=rep("grey",nrow(a))
a$col[which(a$padj<0.05 & a$log2FoldChange>0)]="darkred"
a$col[which(a$padj<0.05 & a$log2FoldChange<0)]="darkblue"

a$hjust <- ifelse(a$log2FoldChange > 0, -0.25, 0.25)
pdf("twister.HNF1A_PE.our_PP_signature.pdf",height=8)
par(mar=c(5.1,8,1,1), xpd=T)
bar=barplot(a$log2FoldChange, xlab="log2FoldChange", col=a$col, cex.lab=1.5, cex.axis=1.5, horiz=T)
text(a$hjust,bar, a$GeneName, cex=1.5)
dev.off()

save.image("DEG.enrichments.Rdata")
