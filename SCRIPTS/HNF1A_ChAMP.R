library("ChAMP")
library("doParallel")
detectCores()

## set analysis directory
dir="/Users/agata/Desktop/StemBANCC/HNF1A/methylation/00.raw_data"
resDir="/Users/agata/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP"

## Sample_sheet.csv has been slightly modified from the original file provided by the core
## to better match the example sample_sheet.csv provided by ChAMP 
## changed column order, lowercase [Data] and [Header], remove [Manifest], uppercase Sentrix_Position to match file names
## Also - renamed the samples, and added Sample_Group as well
## I keep here only the HNF1A project data; HNF4A is moved to a different folder

## load and filter data (done in one step)
myLoad <- champ.load(dir, arraytype="EPIC")

## 6 filters:
## 1) by probe detection value (p>0.01)
## 2) probes with <3 beads in at least 5% of samples per probe
## 3) non-CpG probes
## 4) SNP-related probes - from Zhou's NAR 2016 paper  - we could skip this step, as all samples from the same individual
## 5) filter all multi-hit probes (Nordlund Genome Biology 2013 paper)
## 6) all chr X and chr Y probes - we could skip this step, as all samples from the same individual
# myLoad <- champ.load(dir, arraytype="EPIC", population=FALSE, filterXY=FALSE)

## check distribution of CpG's:
#CpG.GUI(CpG=rownames(myLoad$beta),arraytype="EPIC")

## QC
champ.QC(resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/")
#QC.GUI(beta=myLoad$beta,arraytype="EPIC")

## normalization:
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=4, #plotBMIQ=TRUE, 
                     resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/")

## SVD; add differentiation as batch
myLoad$pd$diff = sapply(strsplit(myLoad$pd$Sample_Name, split="_"), function(x) x[2])
champ.SVD(beta=myNorm,pd=myLoad$pd,resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/")

## run ComBat batch adjustment for the effects of differentiation:
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("diff"))
champ.SVD(beta=myCombat,pd=myLoad$pd,resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/")

## find Differential Methylation Probes
myDMP <- champ.DMP(beta = myCombat,pheno=myLoad$pd$Sample_Group)
## this runs all pairwise contrasts for the 6 levels of Sample_Group; we don't care about some of them

# -----------------------------
#   Start to Compare : MODY_DE, Corr_DE
# Contrast Matrix
# Contrasts
# Levels     pMODY_DE-pCorr_DE
# pCorr_DE                -1
# pMODY_DE                 1
# You have found 15502 significant MVPs with a BH adjusted P-value below 0.05.
#  
# -----------------------------
#   Start to Compare : MODY_PE, Corr_PE
# Contrast Matrix
# Contrasts
# Levels     pMODY_PE-pCorr_PE
# pCorr_PE                -1
# pMODY_PE                 1
# You have found 4801 significant MVPs with a BH adjusted P-value below 0.05.
# 
# -----------------------------
#   Start to Compare : MODY_EN, Corr_EN
# Contrast Matrix
# Contrasts
# Levels     pMODY_EN-pCorr_EN
# pCorr_EN                -1
# pMODY_EN                 1
# You have found 1855 significant MVPs with a BH adjusted P-value below 0.05.

## look at top results for each of the interesting comparisons
head(myDMP[[1]])
head(myDMP[[10]])
head(myDMP[[15]])

## save DMP results
write.table(myDMP[[1]], file=paste0(resDir,"/DE.DMProbes.ChAMP.txt"),sep="\t", quote=F)
write.table(myDMP[[10]], file=paste0(resDir,"/PE.DMProbes.ChAMP.txt"),sep="\t", quote=F)
write.table(myDMP[[15]], file=paste0(resDir,"/EN.DMProbes.ChAMP.txt"),sep="\t", quote=F)

## explore DMRs in the Shiny interface:
#DMP.GUI(DMP=myDMP[[15]],beta=myCombat,pheno=myLoad$pd$Sample_Group)

## Differential Methylation Regions (DMRs)
## extended segments of the genome that show a quantitative alteration in DNA methylation levels between two groups
## 4 methods: Bumphunter, ProbeLasso and DMRcate
## DMRcate is the newest

## champ.DMR() only takes categorical covariates with exactly TWO PHENOTYPES; 
## need to manually select two phenotypes at a time, and input corresponding beta matrix and pheno information
de=grepl("_DE",myLoad$pd$Sample_Group)
pe=grepl("_PE",myLoad$pd$Sample_Group)
en=grepl("_EN",myLoad$pd$Sample_Group)

myDMR_DE <- champ.DMR(beta=myCombat[,de],pheno=myLoad$pd$Sample_Group[de],method="DMRcate" )
#DMRcate detected 484 DMRs with mafcut as= 0.05.
head(myDMR_DE$DMRcateDMR)
myDMR_PE <- champ.DMR(beta=myCombat[,pe],pheno=myLoad$pd$Sample_Group[pe],method="DMRcate" )
# DMRcate detected 65 DMRs with mafcut as= 0.05.
head(myDMR_PE$DMRcateDMR)

myDMR_EN <- champ.DMR(beta=myCombat[,en],pheno=myLoad$pd$Sample_Group[en],method="DMRcate")
head(myDMR_EN$DMRcateDMR)
# DMRcate detected 99 DMRs with mafcut as= 0.05.

## save DMR results:
write.table(myDMR_DE$DMRcateDMR, file=paste0(resDir,"/DE.DMRegions.ChAMP.txt"),sep="\t", quote=F)
write.table(myDMR_PE$DMRcateDMR, file=paste0(resDir,"/PE.DMRegions.ChAMP.txt"),sep="\t", quote=F)
write.table(myDMR_EN$DMRcateDMR, file=paste0(resDir,"/EN.DMRegions.ChAMP.txt"),sep="\t", quote=F)


#DMR.GUI(DMR=myDMR)  ## this GUI doesn't work!

### find differentially methylated blocks (DMBs)
myBlock_DE <- champ.Block(beta=myCombat[,de],pheno=myLoad$pd$Sample_Group[de],arraytype="EPIC")
head(myBlock_DE$Block)
# [bumphunterEngine] Found 1195 bumps.

myBlock_PE <- champ.Block(beta=myCombat[,pe],pheno=myLoad$pd$Sample_Group[pe],arraytype="EPIC")
head(myBlock_PE$Block)
# [bumphunterEngine] Found 958 bumps.

myBlock_EN <- champ.Block(beta=myCombat[,en],pheno=myLoad$pd$Sample_Group[en],arraytype="EPIC")
head(myBlock_EN$Block)
# [bumphunterEngine] Found 1371 bumps.

### save DMBs results
write.table(myBlock_DE$Block, file=paste0(resDir,"/DE.DMBlocks.ChAMP.txt"),sep="\t", quote=F)
write.table(myBlock_PE$Block, file=paste0(resDir,"/PE.DMBlocks.ChAMP.txt"),sep="\t", quote=F)
write.table(myBlock_EN$Block, file=paste0(resDir,"/EN.DMBlocks.ChAMP.txt"),sep="\t", quote=F)


### run GSEA on DMPs and DMRs
myGSEA <- champ.GSEA(beta=myCombat,DMP=myDMP[[1]], DMR=myDMR_DE, arraytype="EPIC",adjPval=0.05, method="gometh")
head(myGSEA$DMP,n=20)
head(myGSEA$DMR,n=20)

myGSEA_PE <- champ.GSEA(beta=myCombat,DMP=myDMP[[10]], DMR=myDMR_PE, arraytype="EPIC",adjPval=0.05, method="gometh")
head(myGSEA_PE$DMP,n=20)
head(myGSEA_PE$DMR,n=20)

myGSEA_EN <- champ.GSEA(beta=myCombat,DMP=myDMP[[15]], DMR=myDMR_EN, arraytype="EPIC",adjPval=0.05, method="gometh")
head(myGSEA_EN$DMP,n=20)
head(myGSEA_EN$DMR,n=20)

### save GSEA results on DMPs and DMRs
write.table(myGSEA$DMP, file=paste0(resDir,"/DE.DMP_GSEA.ChAMP.txt"),sep="\t", quote=F)
write.table(myGSEA$DMR, file=paste0(resDir,"/DE.DMR_GSEA.ChAMP.txt"),sep="\t", quote=F)
write.table(myGSEA_PE$DMP, file=paste0(resDir,"/PE.DMP_GSEA.ChAMP.txt"),sep="\t", quote=F)
write.table(myGSEA_PE$DMR, file=paste0(resDir,"/PE.DMR_GSEA.ChAMP.txt"),sep="\t", quote=F)
write.table(myGSEA_EN$DMP, file=paste0(resDir,"/EN.DMP_GSEA.ChAMP.txt"),sep="\t", quote=F)
write.table(myGSEA_EN$DMR, file=paste0(resDir,"/EN.DMR_GSEA.ChAMP.txt"),sep="\t", quote=F)


### empirical Bayes GSEA:
myebayGSEA_DE <- champ.ebGSEA(beta=myCombat[,de],pheno=myLoad$pd$Sample_Group[de],arraytype="EPIC")
myebayGSEA_PE <- champ.ebGSEA(beta=myCombat[,pe],pheno=myLoad$pd$Sample_Group[pe],arraytype="EPIC")
myebayGSEA_EN <- champ.ebGSEA(beta=myCombat[,en],pheno=myLoad$pd$Sample_Group[en],arraytype="EPIC")

## save empirical Bayes GSEA results:
write.table(myebayGSEA_DE$`Rank(P)` , file=paste0(resDir,"/DE.ebayGSEA_rankP.ChAMP.txt"),sep="\t", quote=F)
write.table(myebayGSEA_DE$`Rank(AUC)`, file=paste0(resDir,"/DE.ebayGSEA_rankAUC.ChAMP.txt"),sep="\t", quote=F)
write.table(myebayGSEA_PE$`Rank(P)`, file=paste0(resDir,"/PE.ebayGSEA_rankP.ChAMP.txt"),sep="\t", quote=F)
write.table(myebayGSEA_PE$`Rank(AUC)`, file=paste0(resDir,"/PE.ebayGSEA_rankAUC.ChAMP.txt"),sep="\t", quote=F)
write.table(myebayGSEA_EN$`Rank(P)`, file=paste0(resDir,"/EN.ebayGSEA_rankP.ChAMP.txt"),sep="\t", quote=F)
write.table(myebayGSEA_EN$`Rank(AUC)`, file=paste0(resDir,"/EN.ebayGSEA_rankAUC.ChAMP.txt"),sep="\t", quote=F)



### find differentially methylated interaction hotspots
myEpiMod_DE <- champ.EpiMod(beta=myCombat[,de],pheno=myLoad$pd$Sample_Group[de],resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/EpiMod_DE")
myEpiMod_PE <- champ.EpiMod(beta=myCombat[,pe],pheno=myLoad$pd$Sample_Group[pe],resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/EpiMod_PE")
myEpiMod_EN <- champ.EpiMod(beta=myCombat[,en],pheno=myLoad$pd$Sample_Group[en],resultsDir="~/Desktop/StemBANCC/HNF1A/methylation/analysis_ChAMP/EpiMod_EN")



