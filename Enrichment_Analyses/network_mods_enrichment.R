
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~                                                                    #~#~#~#~#~
#~#~#~#~#~           Over Representation Analysis                                                      #~#~#~#~#~
#~#~#~#~#~                                                                    #~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


options(stringsAsFactors = F)

# set directory to location of geneInfo files output from the WGCNA code
# and the cell type marker files
# the geneInfo file should have a column for genes and another column for the module assignment 
# these should be called:
# "geneSymbol"
# "moduleColor"

directory="~/"
setwd(directory)

#load libraries
library(WGCNA)

# read in geneInfo files to be analyzed
female_networks=list(AMYdata=read.csv(paste0(directory, "geneInfo_F_amy.csv")),
                     PFCdata=read.csv(paste0(directory, "geneInfo_F_pfc.csv")),
                     HYPdata=read.csv(paste0(directory, "geneInfo_F_hyp.csv")),
                     BLOODdata=read.csv(paste0(directory, "geneInfo_F_bld.csv")))


lapply(female_networks, names)

# read in a list of background genes
# be sure the column containing the genes is named "geneSymbol"
background_genes <- lapply(female_networks, function(x) as.data.frame(na.omit(unique(x[,"geneSymbol"]))))
lapply(background_genes, head)
lapply(background_genes, dim)

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#       perform enrichment analyses with 
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

background_list<-background_genes

data_lists<-female_networks

for (i in 1:length(data_lists)) {
  
  geneR = toupper(data_lists[[i]]$geneSymbol)
  labelR=data_lists[[i]]$moduleColor
  
  # Now run the function!
  testResults = userListEnrichment(geneR, labelR,
                                   fnIn=c(paste0(directory, "gene_lists_for_enrichment_analysis.csv")),
                                   nameOut = paste0("CIE_F_userListEnrichment_", names(data_lists)[i], ".csv"),
                                   useBrainLists=F,
                                   useImmunePathwayLists=F,
                                   useBrainRegionMarkers=F,
                                   usePalazzoloWang=F,
                                   #omitCategories ="grey",
                                   outputCorrectedPvalues=T,
                                   outputGenes=T)
}




# MALES
# read in geneInfo files to be analyzed
male_networks=list(AMYdata=read.csv(paste0(directory, "geneInfo_M_amy.csv")),
                   PFCdata=read.csv(paste0(directory, "geneInfo_M_pfc.csv")),
                   HYPdata=read.csv(paste0(directory, "geneInfo_M_hyp.csv")),
                   BLOODdata=read.csv(paste0(directory, "geneInfo_M_bld.csv")))


lapply(male_networks, names)

# read in a list of background genes
# be sure the column containing the genes is named "geneSymbol"
background_genes <- lapply(male_networks, function(x) as.data.frame(na.omit(unique(x[,"geneSymbol"]))))
lapply(background_genes, head)
lapply(background_genes, dim)

#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#
#       perform enrichment analyses with 
#~#~#~#~#~#~#~#~~#~##~#~#~#~#~#~#~#~~#~#

background_list<-background_genes

data_lists<-male_networks

for (i in 1:length(data_lists)) {
  
  geneR = toupper(data_lists[[i]]$geneSymbol)
  labelR=data_lists[[i]]$moduleColor
  
  # Now run the function!
  testResults = userListEnrichment(geneR, labelR,
                                   fnIn=c(paste0(directory, "gene_lists_for_enrichment_analysis.csv")),
                                   nameOut = paste0("CIE_M_userListEnrichment_", names(data_lists)[i], ".csv"),
                                   useBrainLists=F,
                                   useImmunePathwayLists=F,
                                   useBrainRegionMarkers=F,
                                   usePalazzoloWang=F,
                                   #omitCategories ="grey",
                                   outputCorrectedPvalues=T,
                                   outputGenes=T)
}
