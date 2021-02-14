library(WGCNA)
library(DESeq2)


#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
setwd("~/Yupeng")

# Read the gene counts table 
data0=read.table("gene_counts_table_WGCNA.txt",header=T,row.names=1,sep="\t")
# Normalization with log2(FPKM+1)
sample_metadata = read.csv(file = "sample_info.csv")
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-181],colData = sample_metadata,design = ~ Zone)
mcols(dataExpr_deseq)$basepairs = data0$geneLength
fpkm_matrix = fpkm(dataExpr_deseq)
datExpr = t(log2(fpkm_matrix+1))
# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
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

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM 
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

#===============================================================================
#
#  Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================

# Merge close modules
MEDissThres=0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");

#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
    modules = c(substring(names(merge$oldMEs)[i], 3));
    genes = names(datExpr)
    inModule = is.finite(match(dynamicColors,modules))
    modGenes = genes[inModule]
    modTOM=TOM[inModule,inModule]
    dimnames(modTOM)=list(modGenes,modGenes)
    cyt = exportNetworkToCytoscape(modTOM,
            edgeFile = paste("orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
            nodeFile = paste("orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
            weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
    modules = c(substring(names(merge$newMEs)[i], 3));
    genes = names(datExpr)
    inModule = is.finite(match(dynamicColors,modules))
    modGenes = genes[inModule]
    modTOM=TOM[inModule,inModule]
    dimnames(modTOM)=list(modGenes,modGenes)
    cyt = exportNetworkToCytoscape(modTOM,
            edgeFile = paste("merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
            nodeFile = paste("merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
            weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

#===============================================================================
#
#  Plot the heatmap of module eigen-genes and samples
#
#===============================================================================

library("pheatmap")

# Heatmap of old module eigen-genes and samples
pdf(file="oldMEs.pdf",heigh=80,width=20)
row.names(merge$oldMEs)=names(data0)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()
# Heatmap of new module eigen-genes and samples
pdf(file="newMEs.pdf",heigh=60,width=20)
row.names(merge$newMEs)=names(data0)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()

#=====================================================================================
#
#  Correlation between gene modules and microbial traits
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# Read microbial data as traits
bac_traits = read.table("./WGCNA_180samples/my_modules/b_order_234.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
rownames(MEs) = paste(substr(rownames(MEs), 1, nchar(rownames(MEs))-1), rep(c("1", "2", "3"), 60), sep = "")
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))
# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#=====================================================================================
#
#  Plot heatmap of module-traits relationship
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("module-traits-bacteria-order.pdf", width = 100, height = 30)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bac_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


