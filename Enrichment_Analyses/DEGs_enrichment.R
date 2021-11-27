#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#~#~#~#~#~                                                                    #~#~#~#~#~
#~#~#~#~#~           Over Representation Analysis                                                      #~#~#~#~#~
#~#~#~#~#~                                                                    #~#~#~#~#~
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


options(stringsAsFactors = F)

# Load libraries
library(org.Mm.eg.db)
library(AnnotationDbi)
library(WGCNA) 
library(tidyr)
library(BRETIGEA)
library(imsig)

# path to where differential expression results are located
directory<-"~/"
setwd(directory)

# read in limma DE results
# these files can be output from the limma voom scripts
# ensure that the p value col is called "p.value"
# and that the fold change col is called "logFC"
# and gene symbols are in a col names "mgi_symbol"

#read in limma results files

brain_df<-list(amy_m=read.csv(paste0(directory, "limma_results_amy_M.csv")),
               amy_f=read.csv(paste0(directory, "limma_results_amy_F.csv")),
               hyp_m=read.csv(paste0(directory, "limma_results_hyp_M.csv")),
               hyp_f=read.csv(paste0(directory, "limma_results_hyp_F.csv")),
               pfc_m=read.csv(paste0(directory, "limma_results_pfc_M.csv")),
               pfc_f=read.csv(paste0(directory, "limma_results_pfc_F.csv")))


# get cell type markers from the BRETIGEA R package 
markers_df_mouse_brain
write.csv(markers_df_mouse_brain, paste0(directory, "BRETIGEA_markers_df_mouse_brain.csv"), row.names = F)


# Run enrichment analysis on DEGs using WGCNAs function userlistenrichment and the marker genes from BRETIGEA
# write results to file

for (i in 1:length(brain_df)) {
  deStatus<-ifelse(test = brain_df[[i]][,"p.value"] < 0.05, yes = ifelse(test = brain_df[[i]][,"logFC"] > 0, yes = "Up", no = "Down"), no = "background") 
  symbol<-brain_df[[i]]$mgi_symbol
  userListEnrichment(geneR = symbol, 
                     labelR = deStatus, 
                     fnIn = paste0(directory, "BRETIGEA_markers_df_mouse_brain.csv"), 
                     outputGenes = T,
                     outputCorrectedPvalues = T,
                     nameOut = paste0("cell_type_enrichment_results_BRETIGEA_corrected_", names(brain_df)[[i]], ".csv"))
  
}




#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
# ImSig cell type enrichment analysis
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# get cell type markers from the imsig R package 
sig
write.csv(sig, paste0(directory, "ImSig_markers_immune_cells.csv"), row.names = F)

#read in limma results files

blood_df<-list(blood_m=read.csv(paste0(directory, "limma_results_bld_M.csv")),
               blood_f=read.csv(paste0(directory, "limma_results_bld_F.csv")))


# Run enrichment analysis on DEGs using WGCNAs function userlistenrichment and the marker genes from BRETIGEA
# write results to file

for (i in 1:length(blood_df)) {
  
  deStatus<-ifelse(test = blood_df[[i]][,14] < 0.05, yes = ifelse(test = blood_df[[i]][,11] > 0, yes = "Up", no = "Down"), no = "background") 
  
  symbol<-toupper(blood_df[[i]]$mgi_symbol)
  
  userListEnrichment(geneR = symbol, 
                     labelR = deStatus, 
                     fnIn = paste0(directory, "ImSig_markers_immune_cells.csv"), 
                     outputGenes = T,
                     outputCorrectedPvalues = T,
                     nameOut = paste0("cell_type_enrichment_results_ImSig_corrected_", names(blood_df)[[i]], ".csv"))
  
}





