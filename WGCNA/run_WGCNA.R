rm(list=ls(all=TRUE)) 
library("RColorBrewer")
library("flashClust")
library(WGCNA);
workingDir = "/Users/manoj/Data/ENCODE3/RNASeq/ReRunFinalWGCNA/";
setwd(workingDir); 
options(stringsAsFactors = FALSE);
allowWGCNAThreads(nThreads = 12)



###     Part 1: Data Input     ###


AllData <- data.frame(read.table(file="AllDCCData_e10p5remapdDCCIndx_TPMs_Mtrx.txt", header=TRUE,sep="\t",row.names=1))
dim(AllData)
GeneTPM <- AllData[,1:144]
dim(GeneTPM)
GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) !all(x==0)),]
dim(GeneTPM)
GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) sum(x>0))/(ncol(GeneTPM)-1) >= 0.9 ,]  #gives 18401 with all DCC data
write.table(GeneTPM,"TPMMtrx_Used_Fltrd.txt",quote = F,sep = "\t")


# 90 percent of 144 is 129.6 -- ie., ~14 samples should have non-zero values. eg., 14 liver samples have non-zero TPMs of genes - ie., liver specific genes.

#GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) sum(x>0))/(ncol(GeneTPM)-1) >= 0.8 ,]  #gives 19797 with all DCC data
#GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) sum(x>0))/(ncol(GeneTPM)-1) >= 0.7 ,]  #gives 20941 with all DCC data  ; SFT-b: 4
#GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) sum(x>0))/(ncol(GeneTPM)-1) >= 0.6 ,]  #gives 21955 with all DCC data ; SFT b 6 for max b.
dim(GeneTPM)

#GeneTPM <- GeneTPM[apply(GeneTPM[,-1], 1, function(x) !all(x>2000)),]
#dim(GeneTPM)


GeneTPM_vals <- data.matrix(GeneTPM)
GeneTPM_vals_log <- log2(GeneTPM_vals + 0.00001)

write.csv(GeneTPM_vals_log, file="TPM_Mtrx_ENCODEAll_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_BicorCmplthclust_No0TPM_Log2Trnsfmd_Manual.csv",quote=F)

# REMEMBER -- add "GeneID" in this file

ExpData = read.csv("TPM_Mtrx_ENCODEAll_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_BicorCmplthclust_No0TPM_Log2Trnsfmd_Manual.csv");
dim(ExpData);
names(ExpData);



datExpr0 = as.data.frame(t(ExpData[, -c(1:1)]));
names(datExpr0) = ExpData$GeneID;
rownames(datExpr0) = names(ExpData)[-c(1:1)];

write.table(datExpr0, file="datExpr0_TPM_Mtrx_ENCODEAll_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_BicorCmplthclust_No0TPM_Log2Trnsfmd_Manual.txt",sep = "\t",quote=F)


gsg = goodSamplesGenes(datExpr0, verbose = 9);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "ward.D2");
sizeGrWindow(12,9)
pdf(file = "SampleClustering_allOnlyDCCSamples_Full.pdf", width = 18, height = 5);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2, ylim=c(0,3000))
abline(h = 500, col = "red");
dev.off()

clust = cutreeStatic(sampleTree)
#clust = cutreeStatic(sampleTree, cutHeight = 500, minSize = 45)

table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)




traitData = read.csv("All_AllSmplsFull_ENCODESamples_FinalCodes_SampleTraitsReOrdrd_Clean.csv");
dim(traitData)
names(traitData)
allTraits = traitData;
allTraits = allTraits[, c(2:4) ];

dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

MH_RNA_Samples = rownames(datExpr);
traitRows = match(MH_RNA_Samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "ward.D2")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
col_palette <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
traitColors = numbers2colors(datTraits, signed = TRUE);

pdf(file = "SampleClustering_withTraitHeatmap_all_OnlyDCCSamples_Full.pdf", width = 22, height = 5);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap - signed")
dev.off()
save(datExpr, datTraits, file = "WGCNA_ENCODE_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlySamples_RNA-dataInput.RData")


###     Part 2: Network Construction     ###

lnames = load(file = "WGCNA_ENCODE_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlySamples_RNA-dataInput.RData");
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 9, networkType = "signed hybrid", corFnc="bicor", corOptions = "use = 'p', maxPOutliers = 0.1")

sft$fitIndices

# Plot the results:
pdf(file = "SFT_ScaleIndpndce_MeanCnctvty_all_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor.pdf", width = 8, height = 6);
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="SFT Model Fit,Network = signed hybrid; CorFnc=Bicor; R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# #Below we study how the biological findings depend on the choice of the power.
# # We use the following power for the power adjacency function.
# beta1=xxxxxx   # use the required value to see the slope. 
# Connectivity=softConnectivity(datExpr,power=beta1)-1
# # Let’s create a scale free topology plot. The black curve corresponds to scale free topology and the red curve corresponds to truncated scale free topology.
# par(mfrow=c(1,1))
# scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=F); 



####### 

softPower = 5;




#adjacency = adjacency(datExpr, power = softPower, type = "signed hybrid", corFnc="bicor", corOptions = "use = 'p', maxPOutliers = 0.1");

#TOM = TOMsimilarity(adjacency);
#dissTOM = 1-TOM

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softPower, networkType = "signed hybrid",
                                  corType="bicor", 
                                  TOMType = "signed", verbose = 4);


# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM),method="average");

#hc <- hclust(dist(df,method = "euclidean"), method="centroid") 

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

minModuleSize =  30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid", deepSplit = 3, minClusterSize = minModuleSize, verbose = 4);



table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



MergedNtwrk_name = "WGCNA_ENCODE_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlySamples";


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "ward.D2");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.15
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "geneDendro-3_all_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_BicorCmplthclust_No0TPM.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = MergedNtwrk_name)

####### 

table(mergedColors)
table(moduleColors)


ModColorTableName = paste(paste("Freq_ModuleColors_SignedHybrid",MergedNtwrk_name,sep="_"),"txt",sep=".")

library(plyr)
Freq_ModuleCols = arrange(count(moduleColors),desc(freq))
write.table(Freq_ModuleCols,file=ModColorTableName,quote=FALSE,sep="\t")
dim(Freq_ModuleCols)


#To output expression of the eigen genes of the new merged modules:
write.table(mergedMEs,file="EigenGenesExprssn_table.txt", quote = F, sep = "\t", row.names=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"))


###     Part 3: Relate to Traits     ###

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf(file = paste(paste("ModuleTraitRltnshp_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnly", MergedNtwrk_name, sep="_"),"pdf", sep="."), wi = 8, he = 14)
par(mar = c(5, 12, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = substr(names(MEs),3,length(names(MEs))),
               #yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships: signed TOM; signed hybrid network"))

dev.off()



RqdTrait_Tiss = as.data.frame(datTraits$Tissue);   # Define Required Trait here (Tissue or DevTime)
names(RqdTrait_Tiss) = "Trait_Tiss"

RqdTrait_DevTime = as.data.frame(datTraits$DevTime);   # Define Required Trait here (Tissue or DevTime)
names(RqdTrait_DevTime) = "Trait_DevTime"




Tiss = as.data.frame(datTraits$Tissue);
names(Tiss) = "Tissue"
MET = orderMEs(cbind(MEs, Tiss))


DevTime = as.data.frame(datTraits$DevTime);
names(DevTime) = "DevTime"
MEM = orderMEs(cbind(MEs, DevTime))




modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance_Tiss = as.data.frame(cor(datExpr, RqdTrait_Tiss, use = "p"));
GSPvalue_Tiss = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Tiss), nSamples));
names(geneTraitSignificance_Tiss) = paste("GS.", names(RqdTrait_Tiss), sep="");
names(GSPvalue_Tiss) = paste("p.GS.", names(RqdTrait_Tiss), sep="");

geneTraitSignificance_DevTime = as.data.frame(cor(datExpr, RqdTrait_DevTime, use = "p"));
GSPvalue_DevTime = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_DevTime), nSamples));
names(geneTraitSignificance_DevTime) = paste("GS.", names(RqdTrait_DevTime), sep="");
names(GSPvalue_DevTime) = paste("p.GS.", names(RqdTrait_DevTime), sep="");


datKME=signedKME(datExpr, MEs)

annot = read.table("ENCODE_RNASeq_GeneAnntn_GencM4Feats.txt",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
write.table(probes, file="Probes_list_OnlyDCCSamples_Full",quote = F, row.names = F)
probes2annot = match(probes, annot$GeneID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


datGS.Traits=data.frame(cor(datExpr,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datOutput=data.frame(ProbeID=names(datExpr),
                     annot[probes2annot,],moduleColors,datKME,datGS.Traits)
# save the results in a comma delimited file
write.table(datOutput,"ENCODEAll_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlySamples_ModuleTrait_Results.csv",row.names=F,sep=",")


# ****************************************************************************************************************
# ****************************************************************************************************************


FName_TissTrait = "signedHybridNtwrk_forTiss";
FName_DevTimeTrait = "signedHybridNtwrk_forDevTime";


geneInfo0_Tiss = data.frame(GeneID = probes, geneChrom = annot$Chrom[probes2annot],
                            geneStart = annot$Start[probes2annot],
                            geneEnd = annot$End[probes2annot],
                            geneSymbol = annot$GeneName[probes2annot],
                            geneType = annot$GeneType[probes2annot],
                            geneStrand = annot$Strand[probes2annot],
                            moduleColor = moduleColors,
                            geneTraitSignificance_Tiss,
                            GSPvalue_Tiss)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Tiss, use = "p")));    #change to required trait
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames_Tiss = names(geneInfo0_Tiss)
  geneInfo0_Tiss = data.frame(geneInfo0_Tiss, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0_Tiss) = c(oldNames_Tiss, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_Tiss = order(geneInfo0_Tiss$moduleColor, -abs(geneInfo0_Tiss$GS.Trait));
geneInfo_Tiss = geneInfo0_Tiss[geneOrder_Tiss, ]

write.table(geneInfo_Tiss, file = paste("GeneInfo_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlyModules_", FName_TissTrait, ".txt", sep=""), row.names = FALSE,  sep="\t", quote=F)


geneInfo0_DevTime = data.frame(GeneID = probes,
                            geneChrom = annot$Chrom[probes2annot],
                            geneStart = annot$Start[probes2annot],
                            geneEnd = annot$End[probes2annot],
                            geneSymbol = annot$GeneName[probes2annot],
                            geneType = annot$GeneType[probes2annot],
                            geneStrand = annot$Strand[probes2annot],
                            moduleColor = moduleColors,
                            geneTraitSignificance_DevTime,
                            GSPvalue_DevTime)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Tiss, use = "p")));    #change to required trait
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames_DevTime = names(geneInfo0_DevTime)
  geneInfo0_DevTime = data.frame(geneInfo0_DevTime, geneModuleMembership[, modOrder[mod]],
                              MMPvalue[, modOrder[mod]]);
  names(geneInfo0_DevTime) = c(oldNames_DevTime, paste("MM.", modNames[modOrder[mod]], sep=""),
                            paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_DevTime = order(geneInfo0_DevTime$moduleColor, -abs(geneInfo0_DevTime$GS.Trait));
geneInfo_DevTime = geneInfo0_DevTime[geneOrder_DevTime, ]

write.table(geneInfo_DevTime, file = paste("GeneInfo_All_OnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnlyModules_", FName_DevTimeTrait, ".txt", sep=""), row.names = FALSE,  sep="\t", quote=F)



### Part 3a: Diagnostics ###

# colorh1 = moduleColors
# y is trait (Tissue or DevTime)



datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
RelateEigenGenes = signif(cor(datME, use="p"), 2)
write.table(RelateEigenGenes, file = paste("RelateEigenGenes_signedhybrid_", ".txt", sep=""), row.names = FALSE,  sep="\t", quote=F)

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="ward.D2" )
#par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes - signed hybrid network")


pdf(file = paste(paste("ModuleTraitRltnshpBtwnMdls__AllOnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnly", MergedNtwrk_name, sep="_"),"pdf", sep="."), wi = 12, he = 12)
par(mar = c(5, 12, 3, 3));
#sizeGrWindow(8,9)
plotMEpairs(datME,main = "Relationship between module eigengenes - signed hybrid network")
dev.off()

#sizeGrWindow(8,9)

#Exprssn Clstrng - DevTime Pos Corr -- Signfcnt Enrchd Modules

textMatrix_Pairs = paste("r:", signif(moduleTraitCor, 2), "; p:(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix_Pairs) = dim(moduleTraitCor)
yLabels = substr(names(MEs),3,length(names(MEs)))
rownames(textMatrix_Pairs) <- yLabels

pdf(file = paste(paste("ExprssnClstrng_allOnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnly"),"pdf", sep="."), wi = 16, he = 8)
#mtext("Gene Expression Clustering within Modules that have |r| > 0.5; p<0.0001", outer = FALSE, cex = 1.5)
par(mfrow=c(3,1), mar = c(7,4,7,7) + 0.2, mai=c(0.62,1.82,0.82,0.22))

# which.module="lightgreen";
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#         clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
#         title=paste("Module color:",which.module," - |r| > 0.5 and p<0.0001 in Tissue and DevTime - ForTissue:",textMatrix_Pairs[7,1],"; ForDevTime:",textMatrix_Pairs[7,2],"\n"))
#which.module="midnightblue";
#plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
#        title=paste("Module color:",which.module," - |r| > 0.5 and p<0.0001 in Tissue and DevTime\n"))


which.module="orangered4";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
write.table(datME$which.module)
which.module="coral2";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="green";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="lightcoral";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="lightcyan1";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="navajowhite2";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to DevTime | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))


which.module="darkgreen";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to Tissue | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="palevioletred3";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to Tissue | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="lavenderblush3";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to Tissue | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="lightsteelblue1";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to Tissue | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))
which.module="lightyellow";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
        clabels=c("HB_1_E10","HB_2_E10","HB_1_E11","HB_2_E11","HB_1_E12","HB_2_E12","HB_1_E13","HB_2_E13","HB_1_E14","HB_2_E14","HB_1_E15","HB_2_E15","HB_1_E16","HB_2_E16","HB_1_P0","HB_2_P0","MB_1_E10","MB_2_E10","MB_1_E11","MB_2_E11","MB_1_E12","MB_2_E12","MB_1_E13","MB_2_E13","MB_1_E14","MB_2_E14","MB_1_E15","MB_2_E15","MB_1_E16","MB_2_E16","MB_1_P0","MB_2_P0","FB_1_E10","FB_2_E10","FB_1_E11","FB_2_E11","FB_1_E12","FB_2_E12","FB_1_E13","FB_2_E13","FB_1_E14","FB_2_E14","FB_1_E15","FB_2_E15","FB_1_E16","FB_2_E16","FB_1_P0","FB_2_P0","CF_1_E10","CF_2_E10","CF_1_E11","CF_2_E11","CF_1_E12","CF_2_E12","CF_1_E13","CF_2_E13","CF_1_E14","CF_2_E14","CF_1_E15","CF_2_E15","HT_1_E10","HT_2_E10","HT_1_E11","HT_2_E11","HT_1_E12","HT_2_E12","HT_1_E13","HT_2_E13","HT_1_E14","HT_2_E14","HT_1_E15","HT_2_E15","HT_1_E16","HT_2_E16","HT_1_P0","HT_2_P0","LM_1_E10","LM_2_E10","LM_1_E11","LM_2_E11","LM_1_E12","LM_2_E12","LM_1_E13","LM_2_E13","LM_1_E14","LM_2_E14","LM_1_E15","LM_2_E15","LV_1_E11","LV_2_E11","LV_1_E12","LV_2_E12","LV_1_E13","LV_2_E13","LV_1_E14","LV_2_E14","LV_1_E15","LV_2_E15","LV_1_E16","LV_2_E16","LV_1_P0","LV_2_P0","NT_1_E11","NT_2_E11","NT_1_E12","NT_2_E12","NT_1_E13","NT_2_E13","NT_1_E14","NT_2_E14","NT_1_E15","NT_2_E15","IT_1_E14","IT_2_E14","IT_1_E15","IT_2_E15","IT_1_E16","IT_2_E16","IT_1_P0","IT_2_P0","KD_1_E14","KD_2_E14","KD_1_E15","KD_2_E15","KD_1_E16","KD_2_E16","KD_1_P0","KD_2_P0","LG_1_E14","LG_2_E14","LG_1_E15","LG_2_E15","LG_1_E16","LG_2_E16","LG_1_P0","LG_2_P0","ST_1_E14","ST_2_E14","ST_1_E15","ST_2_E15","ST_1_E16","ST_2_E16","ST_1_P0","ST_2_P0"),
        title=paste("Module color:",which.module," - Correlated to Tissue | Cor&PVal ForTissue:",textMatrix_Pairs[which.module,1],"| Cor&PVal ForDevTime",textMatrix_Pairs[which.module,2]," | nGenes in this module = ",Freq_ModuleCols[Freq_ModuleCols$x==which.module,2],"\n"))

dev.off()


sizeGrWindow(8,7);
which.module="orangered4"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")



module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
pdf(file = paste(paste("ExprssnClstrng_AllMdls_allOnlyDCCSamples_Full_AutoMdlDtct_SgndHybNtwrk_SgndTOM_Bicor_GenesGr0TPMOnly"),"pdf", sep="."), wi = 16, he = 8)
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
dev.off()

# FOR COMBINED MODULE HEATMAP USE: http://pklab.med.harvard.edu/scw2014/WGCNA.html

# 
# #sizeGrWindow(8,9)
# par(mfrow=c(4,1), mar=c(3, 2, 6, 2))
# which.module="saddlebrown";
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#         clabels=c("CF_1_E11","CF_2_E11","FB_1_E11","FB_2_E11","HT_1_E11","HT_2_E11","HB_1_E11","HB_2_E11","LM_1_E11","LM_2_E11","LV_1_E11","LV_2_E11","MB_1_E11","MB_2_E11","NT_1_E11","NT_2_E11","CF_1_E12","CF_2_E12","FB_1_E12","FB_2_E12","HT_1_E12","HT_2_E12","HB_1_E12","HB_2_E12","LM_1_E12","LM_2_E12","LV_1_E12","LV_2_E12","MB_1_E12","MB_2_E12","NT_1_E12","NT_2_E12","CF_1_E13","CF_2_E13","FB_1_E13","FB_2_E13","HT_1_E13","HT_2_E13","HB_1_E13","HB_2_E13","LM_1_E13","LM_2_E13","LV_1_E13","LV_2_E13","MB_1_E13","MB_2_E13","NT_1_E13","NT_2_E13","CF_1_E14","CF_2_E14","FB_1_E14","FB_2_E14","HT_1_E14","HT_2_E14","HB_1_E14","HB_2_E14","IT_1_E14","IT_2_E14","KD_1_E14","KD_2_E14","LM_1_E14","LM_2_E14","LV_1_E14","LV_2_E14","LG_1_E14","LG_2_E14","MB_1_E14","MB_2_E14","NT_1_E14","NT_2_E14","ST_1_E14","ST_2_E14","CF_1_E15","CF_2_E15","FB_1_E15","FB_2_E15","HT_1_E15","HT_2_E15","HB_1_E15","HB_2_E15","IT_1_E15","IT_2_E15","KD_1_E15","KD_2_E15","LM_1_E15","LM_2_E15","LV_1_E15","LV_2_E15","LG_1_E15","LG_2_E15","MB_1_E15","MB_2_E15","NT_1_E15","NT_2_E15","ST_1_E15","ST_2_E15","FB_1_E16","FB_2_E16","HT_1_E16","HT_2_E16","HB_1_E16","HB_2_E16","IT_1_E16","IT_2_E16","KD_1_E16","KD_2_E16","LV_1_E16","LV_2_E16","LG_1_E16","LG_2_E16","MB_1_E16","MB_2_E16","ST_1_E16","ST_2_E16","FB_1_P0","FB_2_P0","HT_1_P0","HT_2_P0","HB_1_P0","HB_2_P0","IT_1_P0","IT_2_P0","KD_1_P0","KD_2_P0","LV_1_P0","LV_2_P0","LG_1_P0","LG_2_P0","MB_1_P0","MB_2_P0","ST_1_P0","ST_2_P0"),rcols=which.module,
#         title=which.module )
# which.module="midnightblue";
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#         clabels=c("CF_1_E11","CF_2_E11","FB_1_E11","FB_2_E11","HT_1_E11","HT_2_E11","HB_1_E11","HB_2_E11","LM_1_E11","LM_2_E11","LV_1_E11","LV_2_E11","MB_1_E11","MB_2_E11","NT_1_E11","NT_2_E11","CF_1_E12","CF_2_E12","FB_1_E12","FB_2_E12","HT_1_E12","HT_2_E12","HB_1_E12","HB_2_E12","LM_1_E12","LM_2_E12","LV_1_E12","LV_2_E12","MB_1_E12","MB_2_E12","NT_1_E12","NT_2_E12","CF_1_E13","CF_2_E13","FB_1_E13","FB_2_E13","HT_1_E13","HT_2_E13","HB_1_E13","HB_2_E13","LM_1_E13","LM_2_E13","LV_1_E13","LV_2_E13","MB_1_E13","MB_2_E13","NT_1_E13","NT_2_E13","CF_1_E14","CF_2_E14","FB_1_E14","FB_2_E14","HT_1_E14","HT_2_E14","HB_1_E14","HB_2_E14","IT_1_E14","IT_2_E14","KD_1_E14","KD_2_E14","LM_1_E14","LM_2_E14","LV_1_E14","LV_2_E14","LG_1_E14","LG_2_E14","MB_1_E14","MB_2_E14","NT_1_E14","NT_2_E14","ST_1_E14","ST_2_E14","CF_1_E15","CF_2_E15","FB_1_E15","FB_2_E15","HT_1_E15","HT_2_E15","HB_1_E15","HB_2_E15","IT_1_E15","IT_2_E15","KD_1_E15","KD_2_E15","LM_1_E15","LM_2_E15","LV_1_E15","LV_2_E15","LG_1_E15","LG_2_E15","MB_1_E15","MB_2_E15","NT_1_E15","NT_2_E15","ST_1_E15","ST_2_E15","FB_1_E16","FB_2_E16","HT_1_E16","HT_2_E16","HB_1_E16","HB_2_E16","IT_1_E16","IT_2_E16","KD_1_E16","KD_2_E16","LV_1_E16","LV_2_E16","LG_1_E16","LG_2_E16","MB_1_E16","MB_2_E16","ST_1_E16","ST_2_E16","FB_1_P0","FB_2_P0","HT_1_P0","HT_2_P0","HB_1_P0","HB_2_P0","IT_1_P0","IT_2_P0","KD_1_P0","KD_2_P0","LV_1_P0","LV_2_P0","LG_1_P0","LG_2_P0","MB_1_P0","MB_2_P0","ST_1_P0","ST_2_P0"),rcols=which.module,
#         title=which.module )
# which.module="darkolivegreen";
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#         clabels=c("CF_1_E11","CF_2_E11","FB_1_E11","FB_2_E11","HT_1_E11","HT_2_E11","HB_1_E11","HB_2_E11","LM_1_E11","LM_2_E11","LV_1_E11","LV_2_E11","MB_1_E11","MB_2_E11","NT_1_E11","NT_2_E11","CF_1_E12","CF_2_E12","FB_1_E12","FB_2_E12","HT_1_E12","HT_2_E12","HB_1_E12","HB_2_E12","LM_1_E12","LM_2_E12","LV_1_E12","LV_2_E12","MB_1_E12","MB_2_E12","NT_1_E12","NT_2_E12","CF_1_E13","CF_2_E13","FB_1_E13","FB_2_E13","HT_1_E13","HT_2_E13","HB_1_E13","HB_2_E13","LM_1_E13","LM_2_E13","LV_1_E13","LV_2_E13","MB_1_E13","MB_2_E13","NT_1_E13","NT_2_E13","CF_1_E14","CF_2_E14","FB_1_E14","FB_2_E14","HT_1_E14","HT_2_E14","HB_1_E14","HB_2_E14","IT_1_E14","IT_2_E14","KD_1_E14","KD_2_E14","LM_1_E14","LM_2_E14","LV_1_E14","LV_2_E14","LG_1_E14","LG_2_E14","MB_1_E14","MB_2_E14","NT_1_E14","NT_2_E14","ST_1_E14","ST_2_E14","CF_1_E15","CF_2_E15","FB_1_E15","FB_2_E15","HT_1_E15","HT_2_E15","HB_1_E15","HB_2_E15","IT_1_E15","IT_2_E15","KD_1_E15","KD_2_E15","LM_1_E15","LM_2_E15","LV_1_E15","LV_2_E15","LG_1_E15","LG_2_E15","MB_1_E15","MB_2_E15","NT_1_E15","NT_2_E15","ST_1_E15","ST_2_E15","FB_1_E16","FB_2_E16","HT_1_E16","HT_2_E16","HB_1_E16","HB_2_E16","IT_1_E16","IT_2_E16","KD_1_E16","KD_2_E16","LV_1_E16","LV_2_E16","LG_1_E16","LG_2_E16","MB_1_E16","MB_2_E16","ST_1_E16","ST_2_E16","FB_1_P0","FB_2_P0","HT_1_P0","HT_2_P0","HB_1_P0","HB_2_P0","IT_1_P0","IT_2_P0","KD_1_P0","KD_2_P0","LV_1_P0","LV_2_P0","LG_1_P0","LG_2_P0","MB_1_P0","MB_2_P0","ST_1_P0","ST_2_P0"),rcols=which.module,
#         title=which.module )
# which.module="darkmagenta";
# plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=dimnames((datExpr[,moduleColors==which.module ]))[[2]] ,rcols=which.module,
#         clabels=c("CF_1_E11","CF_2_E11","FB_1_E11","FB_2_E11","HT_1_E11","HT_2_E11","HB_1_E11","HB_2_E11","LM_1_E11","LM_2_E11","LV_1_E11","LV_2_E11","MB_1_E11","MB_2_E11","NT_1_E11","NT_2_E11","CF_1_E12","CF_2_E12","FB_1_E12","FB_2_E12","HT_1_E12","HT_2_E12","HB_1_E12","HB_2_E12","LM_1_E12","LM_2_E12","LV_1_E12","LV_2_E12","MB_1_E12","MB_2_E12","NT_1_E12","NT_2_E12","CF_1_E13","CF_2_E13","FB_1_E13","FB_2_E13","HT_1_E13","HT_2_E13","HB_1_E13","HB_2_E13","LM_1_E13","LM_2_E13","LV_1_E13","LV_2_E13","MB_1_E13","MB_2_E13","NT_1_E13","NT_2_E13","CF_1_E14","CF_2_E14","FB_1_E14","FB_2_E14","HT_1_E14","HT_2_E14","HB_1_E14","HB_2_E14","IT_1_E14","IT_2_E14","KD_1_E14","KD_2_E14","LM_1_E14","LM_2_E14","LV_1_E14","LV_2_E14","LG_1_E14","LG_2_E14","MB_1_E14","MB_2_E14","NT_1_E14","NT_2_E14","ST_1_E14","ST_2_E14","CF_1_E15","CF_2_E15","FB_1_E15","FB_2_E15","HT_1_E15","HT_2_E15","HB_1_E15","HB_2_E15","IT_1_E15","IT_2_E15","KD_1_E15","KD_2_E15","LM_1_E15","LM_2_E15","LV_1_E15","LV_2_E15","LG_1_E15","LG_2_E15","MB_1_E15","MB_2_E15","NT_1_E15","NT_2_E15","ST_1_E15","ST_2_E15","FB_1_E16","FB_2_E16","HT_1_E16","HT_2_E16","HB_1_E16","HB_2_E16","IT_1_E16","IT_2_E16","KD_1_E16","KD_2_E16","LV_1_E16","LV_2_E16","LG_1_E16","LG_2_E16","MB_1_E16","MB_2_E16","ST_1_E16","ST_2_E16","FB_1_P0","FB_2_P0","HT_1_P0","HT_2_P0","HB_1_P0","HB_2_P0","IT_1_P0","IT_2_P0","KD_1_P0","KD_2_P0","LV_1_P0","LV_2_P0","LG_1_P0","LG_2_P0","MB_1_P0","MB_2_P0","ST_1_P0","ST_2_P0"),rcols=which.module,
#         title=which.module )






signif(cor(DevTime,datME, use="p"),2)
cor.test(DevTime, datME$MEturquoise)
p.values = corPvalueStudent(cor(DevTime,datME, use="p"), nSamples = length(DevTime))

GS1=as.numeric(cor(DevTime,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)
ModuleSignificance

#sizeGrWindow(8,7)
pdf(file = paste(paste("ModuleSgnfcnce_AvgGeneSgnfcne_allOnlyDCCSamples_Full"),"pdf", sep="."), wi = 12, he = 12)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors)
dev.off()


###     Part 4: Annotation     ###

# Do this separately... 






###     Part 5: Visualization     ###


# This takes a lot of time, and gives an error: 
# plotTOM = dissTOM^7;
# diag(plotTOM) = NA;
# sizeGrWindow(9,9)
#  TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 1000
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "ward.D2")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
#TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")



pdf(file = paste(paste("Heatmap_TraitTissue_Module_Clstrng_allOnlyDCCSamples_Full"),"pdf", sep="."), wi = 12, he = 12)
#sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),plotHeatmaps = FALSE)

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)


pdf(file = paste(paste("Heatmap_TraitDevTime_Module_Clstrng_allOnlyDCCSamples_Full"),"pdf", sep="."), wi = 12, he = 12)
#sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEM, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEM, "Eigengene dendrogram", marDendro = c(0,4,2,0),plotHeatmaps = FALSE)

par(cex = 1.0)
plotEigengeneNetworks(MEM, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)





###     Part 6: Export Network     ###


####### MH STOP HERE ########



## FOR VisANT
# 
# module_rqd = "magenta";
# probes = names(datExpr)
# inModule = (moduleColors==module_rqd);
# modProbes = probes[inModule];
# # Select the corresponding Topological Overlap
# modTOM = TOM[inModule, inModule];
# dimnames(modTOM) = list(modProbes, modProbes)
# # Export the network into an edge list file VisANT can read
# vis = exportNetworkToVisANT(modTOM,
#                             file = paste("VisANTInput-", module_rqd, ".txt", sep=""),
#                             weighted = TRUE,
#                             threshold = 0,
#                             probeToGene = data.frame(annot$GeneID, annot$GeneName) )
# nTop = 30;
# IMConn = softConnectivity(datExpr[, modProbes]);
# top = (rank(-IMConn) <= nTop)
# vis = exportNetworkToVisANT(modTOM[top, top],
#                             file = paste("VisANTInput-", module_rqd, "-top30.txt", sep=""),
#                             weighted = TRUE,
#                             threshold = 0,
#                             probeToGene = data.frame(annot$GeneID, annot$GeneName) )
# 
# 
# 
# 
# ## FOR CYTOSCAPE 
# 
# # Select modules
# modules = c("magenta", "brown");
# # Select module probes
# probes = names(datExpr)
# inModule = is.finite(match(moduleColors, modules));
# modProbes = probes[inModule];
# modGenes = annot$GeneName[match(modProbes, annot$GeneID)];
# # Select the corresponding Topological Overlap
# modTOM = TOM[inModule, inModule];
# dimnames(modTOM) = list(modProbes, modProbes)
# # Export the network into edge and node list files Cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
#                                weighted = TRUE, threshold = 0.02, nodeNames = modProbes, altNodeNames = modGenes, nodeAttr = moduleColors[inModule]);
# 
# 
