
#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#
#
#     WGCNA
#
#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#

options(stringsAsFactors = F)

library(WGCNA)
library(flashClust)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(DescTools)

# path to where normalized gene expression data are located
directory<-"~/"
setwd(directory)

#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#
#
#   AMY F
#
#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#

# Read in voom transformed expression dataset
# this can be output from the DEG scripts (write the output of the voom command to file)

A2Data <- read.delim(paste0(directory, "voom_amy_F.txt"))
dim(A2Data)
datExpr=as.data.frame(t(A2Data))
dim(datExpr)
head(names(datExpr))
#[1] "ENSMUSG00000000049.11" "ENSMUSG00000000056.7"  "ENSMUSG00000000058.6" 
#[4] "ENSMUSG00000000078.7"  "ENSMUSG00000000085.16" "ENSMUSG00000000088.7"
head(row.names(datExpr))
rownames(datExpr)=substr(rownames(datExpr), 12,21)
head(row.names(datExpr))
#[1] "21.F.V.AMY" "22.F.V.AMY" "23.F.V.AMY" "24.F.V.AMY" "25.F.V.AMY" "26.F.C.AMY"

# Set soft-thresholding powers from 1 to 20.
powers = c(1:20)
# Network topology analysis
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

# Select an appropriate soft-thresholding power.
softPower = sft$powerEstimate
#9
# Calculate the adjacencies.
adjacency = adjacency(datExpr, power = softPower, type='signed')
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed")
# Calculate the corresponding dissimilarity.
dissTOM <- 1-TOM
# Call the hierarchical clustering function.
geneTree = flashClust(as.dist(dissTOM), method = "average")
# Set the minimum module size.
minModuleSize = 100
# Set the cutting height.
detectCutHeight = 0.99
# Module identification using dynamic tree cut.
dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=detectCutHeight, deepSplit = T, minClusterSize = minModuleSize);
# Display module size for each module.
table(dynamicMods)
# Convert numeric labels into colors.
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath.
plotDendroAndColors(dendro=geneTree, colors=dynamicColors, groupLabels="35.99T", rowText=dynamicColors, cex.rowText = 0.5, dendroLabels = F, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = paste0("Gene dendrogram and module colors"))

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#  Find eigengenes and display relationship between them 
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# # Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#  Merge close modules
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

MEDissThres = 0.25
# Plot the cut line into the dendrogram
#abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#sizeGrWindow(12, 9)
pdf(file = paste0(directory, "network_dendrogram_F_amy_merged_mods.pdf"), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# # Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# # Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# # Plot the result
file_name<-paste0(directory,"/ME_dendrogram_F_amy_merged_mods.pdf")
pdf(file=file_name, width = 10, height = 7)# Set graphical parameters.
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#  Correlate merged eigengenes with traits 
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs = mergedMEs

# read in trait data
allTraits=read.csv("~/Downloads/datTraits.csv")
#be sure we have two digits for mouse number
allTraits$Mouse..=Format(x = allTraits$Mouse.., digits = 0, leading = "00")
samples = rownames(datExpr);
traitRows = match(substr(x = samples, start = 0, stop = 2), allTraits$Mouse..);
datTraits = allTraits[traitRows, , drop=F];
rownames(datTraits) = allTraits[traitRows, 1];

# correlations with EtOH intake (g/kg)

moduleTraitCorCIE_baseline = cor(MEs, datTraits$BASELINEaverage, use = "p");
moduleTraitPvalueCIE_baseline = corPvalueStudent(moduleTraitCorCIE_baseline, nSamples);
moduleTraitPvalue_adjustedCIE_baseline = p.adjust(moduleTraitPvalueCIE_baseline, method = "fdr", n = length(moduleTraitPvalueCIE_baseline))
moduleTraitCorCIE_baseline[which (moduleTraitPvalue_adjustedCIE_baseline<.05)]
moduleTraitCorCIE_baseline[which (moduleTraitPvalueCIE_baseline<.05), row.names=T]
#no modules correlated with baseline drinking levels
CIE_baseline = paste(signif(moduleTraitCorCIE_baseline, 2), "\n(",
                    signif(moduleTraitPvalueCIE_baseline, 1), ")", sep = "")

moduleTraitCorCIE1 = cor(MEs, datTraits$CIE1avergae, use = "p");
moduleTraitPvalueCIE1 = corPvalueStudent(moduleTraitCorCIE1, nSamples);
moduleTraitPvalue_adjustedCIE1 = p.adjust(moduleTraitPvalueCIE1, method = "fdr", n = length(moduleTraitPvalueCIE1))
moduleTraitCorCIE1[which (moduleTraitPvalue_adjustedCIE1<.05)]
moduleTraitCorCIE1[which (moduleTraitPvalueCIE1<.05), row.names=T]
#no modules correlated with first CIE drinking
CIE1 = paste(signif(moduleTraitCorCIE1, 2), "\n(",
                    signif(moduleTraitPvalueCIE1, 1), ")", sep = "")


moduleTraitCorCIE2 = cor(MEs, datTraits$CIE2avergae, use = "p");
moduleTraitPvalueCIE2 = corPvalueStudent(moduleTraitCorCIE2, nSamples);
moduleTraitPvalue_adjustedCIE2 = p.adjust(moduleTraitPvalueCIE2, method = "fdr", n = length(moduleTraitPvalueCIE2))
moduleTraitCorCIE2[which (moduleTraitPvalue_adjustedCIE2<.05)]
moduleTraitCorCIE2[which (moduleTraitPvalueCIE2<.05), row.names=T]
#MEcyan    MEbrown 
#-0.5156754  0.4807770
CIE2 = paste(signif(moduleTraitCorCIE2, 2), "\n(",
                    signif(moduleTraitPvalueCIE2, 1), ")", sep = "")


moduleTraitCorCIE3 = cor(MEs, datTraits$CIE3avergae, use = "p");
moduleTraitPvalueCIE3 = corPvalueStudent(moduleTraitCorCIE3, nSamples);
moduleTraitPvalue_adjustedCIE3 = p.adjust(moduleTraitPvalueCIE3, method = "fdr", n = length(moduleTraitPvalueCIE3))
moduleTraitCorCIE3[which (moduleTraitPvalue_adjustedCIE3<.05)]
moduleTraitCorCIE3[which (moduleTraitPvalueCIE3<.05), row.names=T]
#no modules correlated with third CIE drinking
CIE3 = paste(signif(moduleTraitCorCIE3, 2), "\n(",
                    signif(moduleTraitPvalueCIE3, 1), ")", sep = "")


moduleTraitCor_lastCIE = cor(MEs, datTraits$CIE4avergae, use = "p");
moduleTraitPvalue_lastCIE = corPvalueStudent(moduleTraitCor_lastCIE, nSamples);
moduleTraitPvalue_adjusted_lastCIE = p.adjust(moduleTraitPvalue_lastCIE, method = "fdr", n = length(moduleTraitPvalue_lastCIE))
moduleTraitCor_lastCIE[which (moduleTraitPvalue_adjusted_lastCIE<.05)]
moduleTraitCor_lastCIE[which (moduleTraitPvalue_lastCIE<.05), row.names=T]
#no modules correlated with lastCIE drinking
lastCIE = paste(signif(moduleTraitCor_lastCIE, 2), "\n(",
                    signif(moduleTraitPvalue_lastCIE, 1), ")", sep = "")


# correlations with EtOH Preference (no units)
moduleTraitCorPref_baseline = cor(MEs, datTraits$averagePrefBaseline, use = "p");
moduleTraitPvaluePref_baseline = corPvalueStudent(moduleTraitCorPref_baseline, nSamples);
moduleTraitPvalue_adjustedPref_baseline = p.adjust(moduleTraitPvaluePref_baseline, method = "fdr", n = length(moduleTraitPvaluePref_baseline))
moduleTraitCorPref_baseline[which (moduleTraitPvalue_adjustedPref_baseline<.05)]
# MEred  MEdarkgreen
#[1] -0.6974097  0.7773462
moduleTraitCorPref_baseline[which (moduleTraitPvaluePref_baseline<.05), row.names=T]
#MEturquoise   MEgreenyellow      MEdarkgrey           MEred MEpaleturquoise          MEblue     MEdarkgreen 
#-0.5506992      -0.5593500      -0.4758715      -0.6974097      -0.4894591       0.5164721       0.7773462 
#MEdarkorange 
#0.4707474 
Pref_baseline = paste(signif(moduleTraitCorPref_baseline, 2), "\n(",
                    signif(moduleTraitPvaluePref_baseline, 1), ")", sep = "")


moduleTraitCorPref1 = cor(MEs, datTraits$averagePrefCIE1, use = "p");
moduleTraitPvaluePref1 = corPvalueStudent(moduleTraitCorPref1, nSamples);
moduleTraitPvalue_adjustedPref1 = p.adjust(moduleTraitPvaluePref1, method = "fdr", n = length(moduleTraitPvaluePref1))
moduleTraitCorPref1[which (moduleTraitPvalue_adjustedPref1<.05)]
moduleTraitCorPref1[which (moduleTraitPvaluePref1<.05), row.names=T]
#no modules correlated with first Pref drinking
Pref1 = paste(signif(moduleTraitCorPref1, 2), "\n(",
                    signif(moduleTraitPvaluePref1, 1), ")", sep = "")


moduleTraitCorPref2 = cor(MEs, datTraits$averagePrefCIE2, use = "p");
moduleTraitPvaluePref2 = corPvalueStudent(moduleTraitCorPref2, nSamples);
moduleTraitPvalue_adjustedPref2 = p.adjust(moduleTraitPvaluePref2, method = "fdr", n = length(moduleTraitPvaluePref2))
moduleTraitCorPref2[which (moduleTraitPvalue_adjustedPref2<.05)]
moduleTraitCorPref2[which (moduleTraitPvaluePref2<.05), row.names=T]
#MEgrey 
#0.518238
Pref2 = paste(signif(moduleTraitCorPref2, 2), "\n(",
                    signif(moduleTraitPvaluePref2, 1), ")", sep = "")

moduleTraitCorPref3 = cor(MEs, datTraits$averagePrefCIE3, use = "p");
moduleTraitPvaluePref3 = corPvalueStudent(moduleTraitCorPref3, nSamples);
moduleTraitPvalue_adjustedPref3 = p.adjust(moduleTraitPvaluePref3, method = "fdr", n = length(moduleTraitPvaluePref3))
moduleTraitCorPref3[which (moduleTraitPvalue_adjustedPref3<.05)]
moduleTraitCorPref3[which (moduleTraitPvaluePref3<.05), row.names=T]
#no modules correlated with third Pref drinking
Pref3 = paste(signif(moduleTraitCorPref3, 2), "\n(",
                    signif(moduleTraitPvaluePref3, 1), ")", sep = "")


moduleTraitCor_lastPref = cor(MEs, datTraits$averagePrefCIE4, use = "p");
moduleTraitPvalue_lastPref = corPvalueStudent(moduleTraitCor_lastPref, nSamples);
moduleTraitPvalue_adjusted_lastPref = p.adjust(moduleTraitPvalue_lastPref, method = "fdr", n = length(moduleTraitPvalue_lastPref))
moduleTraitCor_lastPref[which (moduleTraitPvalue_adjusted_lastPref<.05)]
moduleTraitCor_lastPref[which (moduleTraitPvalue_lastPref<.05), row.names=T]
#no modules correlated with lastPref drinking
lastPref = paste(signif(moduleTraitCor_lastPref, 2), "\n(",
                    signif(moduleTraitPvalue_lastPref, 1), ")", sep = "")


# correlations with weight (g)

moduleTraitCorfirstWeight = cor(MEs, datTraits$WeightStart, use = "p");
moduleTraitPvaluefirstWeight = corPvalueStudent(moduleTraitCorfirstWeight, nSamples);
moduleTraitPvalue_adjustedfirstWeight = p.adjust(moduleTraitPvaluefirstWeight, method = "fdr", n = length(moduleTraitPvaluefirstWeight))
moduleTraitCorfirstWeight[which (moduleTraitPvalue_adjustedfirstWeight<.05)]
moduleTraitCorfirstWeight[which (moduleTraitPvaluefirstWeight<.05), row.names=T]
#MEorange 
#-0.5094924 
firstWeight = paste(signif(moduleTraitCorfirstWeight, 2), "\n(",
                     signif(moduleTraitPvaluefirstWeight, 1), ")", sep = "")


moduleTraitCorlastWeight = cor(MEs, datTraits$WeightEnd, use = "p");
moduleTraitPvaluelastWeight = corPvalueStudent(moduleTraitCorlastWeight, nSamples);
moduleTraitPvalue_adjustedlastWeight = p.adjust(moduleTraitPvaluelastWeight, method = "fdr", n = length(moduleTraitPvaluelastWeight))
moduleTraitCorlastWeight[which (moduleTraitPvalue_adjustedlastWeight<.05)]
moduleTraitCorlastWeight[which (moduleTraitPvaluelastWeight<.05), row.names=T]
#MEgreenyellow 
#0.4711776 
lastWeight = paste(signif(moduleTraitCorlastWeight, 2), "\n(",
                    signif(moduleTraitPvaluelastWeight, 1), ")", sep = "")

# Correlations with dependent status

datTraits$Group<-as.numeric(as.factor(datTraits$Group))
moduleTraitCor_Group = cor(MEs, datTraits$Group, use = "p");
moduleTraitPvalue_Group = corPvalueStudent(moduleTraitCor_Group, nSamples);
moduleTraitPvalue_adjusted_Group = p.adjust(moduleTraitPvalue_Group, method = "fdr", n = length(moduleTraitPvalue_Group))
moduleTraitCor_Group[which (moduleTraitPvalue_adjusted_Group<.05)]
moduleTraitCor_Group[which (moduleTraitPvalue_Group<.05), row.names=T]
#no modules correlated with Dependence status
Dependence = paste(signif(moduleTraitCor_Group, 2), "\n(",
                   signif(moduleTraitPvalue_Group, 1), ")", sep = "")


# Will display correlations and their p-values
dataMatrix2 = cbind(moduleTraitCorCIE_baseline, moduleTraitCorCIE1, moduleTraitCorCIE2, moduleTraitCorCIE3, moduleTraitCor_lastCIE, moduleTraitCorPref_baseline, moduleTraitCorPref1, moduleTraitCorPref2, moduleTraitCorPref3, moduleTraitCor_lastPref, moduleTraitCorfirstWeight, moduleTraitCorlastWeight, moduleTraitCor_Group)
textMatrix2 = cbind(CIE_baseline, CIE1, CIE2, CIE3, lastCIE, Pref_baseline, Pref1, Pref2, Pref3, lastPref, firstWeight, lastWeight, Dependence)
dim(textMatrix2) = dim(moduleTraitCor_Dependence)

p_matrix = cbind(moduleTraitPvalueCIE_baseline, moduleTraitPvalueCIE1, moduleTraitPvalueCIE2, moduleTraitPvalueCIE3, moduleTraitPvalue_lastCIE, moduleTraitPvaluePref_baseline, moduleTraitPvaluePref1, moduleTraitPvaluePref2, moduleTraitPvaluePref3, moduleTraitPvalue_lastPref, moduleTraitPvaluefirstWeight, moduleTraitPvaluelastWeight, moduleTraitPvalue_Group)
colnames(p_matrix)=colnames(textMatrix2)
p_matrix = -log(x = p_matrix, base = 10)
colnames(dataMatrix2)=colnames(textMatrix2)
#write.csv(p_matrix, paste0(directory,"CIE_F_amy_module_trait_P_table.csv"))
#write.csv(dataMatrix2, paste0(directory,"CIE_F_amy_module_trait_corr_table.csv"))


file_name<-paste0(directory,"ME_trait_correlations_amy_F_merged_MEs.pdf")
pdf(file = file_name, width = 12, height = 25)
par(mar=c(4,9,1.1,2.1)) #set margins: bottom, left, top and right
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = dataMatrix2,
               xLabels = colnames(textMatrix2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste0("Module - Trait relationship"))


dev.off()


#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#       
#       get transcript annotation
#
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
library(biomaRt)
library(dplyr)
ensembl<- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl",host="www.ensembl.org") #you can use this if using R on your local machine
#use this is using POD
#ensembl<- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
annot<-getBM(c("ensembl_gene_id_version","ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype", "description"), values= names(datExpr), mart=ensembl)
idx <- match(gsub("(.*)(\\..*)","\\1",names(datExpr)),annot$ensembl_gene_id)
res <- data.frame(cbind(names(datExpr), annot[idx,]))
dim(res)
names(res)
annot<-res

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#
#       write gene info file 
#
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~~#~#
#
#      gene significance based on correlation to trait and module membership
#
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~~#~#

MEs = mergedMEs
Name = row.names(datExpr)
Treatment= substr(Name, 6,6)
allTraits=as.data.frame(cbind(Name, Treatment), stringsAsFactors = F)
# Form a data frame analogous to expression data that will hold the clinical traits.
samples = rownames(datExpr);
traitRows = match(samples, allTraits$Name);
datTraits = allTraits[traitRows, -1, drop=F];
rownames(datTraits) = allTraits[traitRows, 1];
# Define variable treatment containing the treatment column of datTrait
Treatment = as.data.frame(as.numeric(as.factor(datTraits$Treatment)))
names(Treatment) = "Treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Treatment, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Treatment), sep="");
names(GSPvalue) = paste("p.GS.", names(Treatment), sep="");

module<-substr(x = row.names(moduleTraitCor_Dependence)[which (moduleTraitPvalue_Dependence<.05)],
               start = 3, stop = nchar(row.names(moduleTraitCor_Dependence)[which (moduleTraitPvalue_Dependence<.05)]))


for ( i in 1:length(module)) {
  column = match(module[[i]], modNames);
  moduleGenes = mergedColors==module[[i]];
  
  pdf(file=paste0(directory, module[[i]], "_corr_vs_MM_F_amy.pdf"))
  
  par(mfrow = c(1,1));
  
  verboseScatterplot(geneModuleMembership[moduleGenes, column],
                     geneTraitSignificance[moduleGenes, 1],
                     xlab = paste("Module Membership in ", module[[i]], " module"),
                     ylab = "Gene correlation with dependence status",
                     main = paste("Module membership vs. gene significance"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module[[i]])
  
  text(geneModuleMembership[moduleGenes, column],
       geneTraitSignificance[moduleGenes, 1], 
       labels = annot[match(x=row.names(geneModuleMembership)[moduleGenes], table = annot$ensembl_gene_id_version),4], pos=4, cex=0.5)
  
  dev.off()
  
  file_name=paste0(directory, module[[i]], "_genes_F_amy.csv")
  
  genes<-annot[match(x=row.names(geneModuleMembership)[moduleGenes], table = annot$ensembl_gene_id_version),c(2,4,10)]
  
  #write.csv(genes, file = file_name)
  
}

probes = gsub("(.*)(\\..*)","\\1",names(datExpr))
probes2annot = match(probes, annot$ensembl_gene_id)

#We now create a data frame holding the following information for all probes: 
#probe ID, gene symbol, module color, gene significance for alc pref, and module membership and p-values in all modules. 
#The modules will be ordered by their significance for alc pref, with the most significant ones to the left.

# Calculate intramodular connectivity.
connectivity=intramodularConnectivity(adjacency, mergedColors)


geneInfo0 = data.frame(ensembl_gene_id_version = probes,
                       geneSymbol = annot$mgi_symbol[probes2annot],
                       geneDescription = annot$description[probes2annot],
                       gene_biotype = annot$gene_biotype[probes2annot],
                       moduleColor = mergedColors,
                       connectivity[,1:2], 
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for treatment
modOrder = order(-abs(cor(MEs, Treatment, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Treatment));
geneInfo = geneInfo0[geneOrder, ]
geneInfo$DE = apply(X = geneInfo[,c(grep(pattern = ".p$", x = names(geneInfo)))], MARGIN = 1, function(x) sum(x<=0.05))

file_name<-paste0(directory, "geneInfo_F_amy.csv")

#write.csv(geneInfo, file = file_name)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#     Network visualization
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Cytoscape is a powerful program for network analysis and visualization. Here we demonstrate how to generate an edge file and a node file for network visualization. Please refer Cytoscape documentation for details on how to input these files into Cytoscape.
# Select a module or modules 
module = c("darkorange",
           "cyan",
           "darkgrey",
           "greenyellow",
           "paleturquoise",
           "red")

# Get all probe IDs.
ensemblIDs = colnames(datExpr)

# Find probes 
inModule = is.finite(match(dynamicColors, module))

# Select probes in the modules
mod_ensemblIDs = ensemblIDs[inModule]

# Read in the annotation file.
ensembl<- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl",host="www.ensembl.org") #you can use this if using R on your local machine
#use this is using POD
#ensembl<- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
annot<-getBM(c("ensembl_gene_id_version","mgi_symbol"), values= ensemblIDs, mart=ensembl)
idx <- match(ensemblIDs,annot$ensembl_gene_id_version)
res <- data.frame(cbind(ensemblIDs, annot[idx,]))
dim(res)
names(res)
annot<-res

# Get the corresponding gene symbols.
modGenes = annot$mgi_symbol[match(mod_ensemblIDs, annot$ensemblIDs)]

# Select the corresponding Topological Overlap.
modTOM = TOM[inModule, inModule]

# Assign row and column names.
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read. 
# Note, the parameter, threshold is adjustable for desired visualization.
cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInput_edges_CIE_F_amy", paste(module, collapse="_"), ".txt", sep=""),nodeFile = paste("CytoscapeInput_nodes_CIE_F_amy", paste(module, collapse="_"), ".txt", sep=""),weighted = TRUE,threshold = 0.02,nodeNames = mod_ensemblIDs, altNodeNames = modGenes, nodeAttr = dynamicColors[inModule])
