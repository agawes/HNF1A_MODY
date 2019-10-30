options(width=200)

setwd("~/Desktop/StemBANCC/HNF1A/RNA-seq bulk/ANALYSIS/")
load("MDS_PCA.Rdata")
load("HNF1A.WGCNA.Rdata")

####### DGCA #####
### probably too small sample size for differential co-expression ####
library(DGCA, quietly = TRUE)
EN_des_DGCA=data.matrix(data.frame(MODY=-as.numeric(design_EN$HNF1A)+2, Corrected=as.numeric(design_EN$HNF1A)-1))
ddcor_res = ddcorAll(inputMat = EN_wgcna_vst, design = EN_des_DGCA,
                     compare = c("MODY","Corrected"),
                     adjust = "none", nPerm = 0, nPairs = 100)
# Error: vector memory exhausted (limit reached?)

## for a single gene:
g="HNF1A"
ens_g=rownames(gene_ann[which(gene_ann$GeneName==g),,drop=F])
ddcor_res = ddcorAll(inputMat = EN_wgcna_vst, design = EN_des_DGCA,
                     compare = c("MODY","Corrected"),
                     adjust = "none", nPerm = 0, splitSet = ens_g)

ddcor_res$GeneName=sapply(ddcor_res$Gene1, function(x) as.character(gene_ann[x,]))
head(ddcor_res)
write.table(ddcor_res,paste0(g,".diff_cor_genes.EN.DGCA.txt"),sep="\t",quote=F)

plotCors(inputMat = EN_wgcna_vst, design = EN_des_DGCA,
         compare = c("MODY","Corrected"), geneA = "ENSG00000254647.2", geneB = "ENSG00000047056.10")

plotCors(inputMat = EN_wgcna_vst, design = EN_des_DGCA,
         compare = c("MODY","Corrected"), geneA = "ENSG00000066032.14", geneB = "ENSG00000135100.13")

