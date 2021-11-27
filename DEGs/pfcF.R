#DE Analysis TagSeq Data (RNA seq) data

# Set options
options(stringsAsFactors = F)

# Load libraries
library(edgeR)
library(limma)
library(biomaRt)
library(dplyr)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Read in count data
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# set directory to location of count files which can be downloaded from GEO (GSE176122)

directory<-"~/Downloads/GSE176122_RAW/"
setwd(directory)

# Read in the count files 

pfc_files<-list.files(directory)[grep(pattern = "F.*PFC|F.*pfc", x = list.files(directory))]

x_pfc <- readDGE(pfc_files, columns= c(1,2)) 
class(x_pfc)
dim(x_pfc)
#[1] 55481    19

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize sample information
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Assign sample-level information related to the experimental design to the columns of the counts matrix [treatment group (CIE vs Air), sex (M and F), tissue (blood, pfcgdala, prefrontal cortex, hypothalamus)] 
samplenames <- substring(colnames(x_pfc), 1, 21)
colnames(x_pfc) <- samplenames
treatment <- as.factor(substr(samplenames, 17, 17))

x_pfc$samples$treatment <- treatment

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize gene annotations
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
geneid <- rownames(x_pfc) 
geneid <- gsub("(.*)(\\..*)", "\\1", geneid)

ensembl<- useEnsembl(biomart="ensembl", 
                     dataset="mmusculus_gene_ensembl",
                     mirror = "useast") 

genes<-getBM(attributes = c("ensembl_gene_id_version","ensembl_gene_id", "mgi_symbol", "external_gene_name", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"),
             filters = "ensembl_gene_id",
             values= geneid, 
             mart=ensembl)

dim(genes)
# [1] 52324     9
genes <- genes[!duplicated(genes$ensembl_gene_id),]
dim(genes)
#[1] 52322     9

#the gene order needs to be the same in both the annotation and the data matrix
genes<-genes[match(x = gsub("(.*)(\\..*)","\\1",rownames(x_pfc)), table = genes$ensembl_gene_id),]
x_pfc$genes<-genes

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      remove ERCC counts from the count files
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

ercc.exprs <- grep(pattern = "GERCC", x = rownames(x_pfc))
# remove ERCCs
x_pfc <-x_pfc[-ercc.exprs,, keep.lib.sizes=FALSE]
dim(x_pfc)
## [1] 55389    19


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      Remove genes that are lowly expressed
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

keep.exprs <- rowSums(cpm(x_pfc)>1)>=8
x_pfc <- x_pfc[keep.exprs,, keep.lib.sizes=FALSE]
dim(x_pfc)
#[1] 15090   19

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Normalize gene expression distributions
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

x_pfc <- calcNormFactors(x_pfc, method = "TMM")

x_pfc$samples$norm.factors
# [1] 0.9630629 0.9861514 0.9507392 0.9967883 1.0202023 1.0240919 0.9874015 0.9690423 0.9705621
# [10] 1.0066512 1.0051369 1.0020730 0.9972018 1.0306320 1.0156726 1.0048801 1.0091615 1.0051529
# [19] 1.0614541

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#         Differential expression analysis
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      Create a design matrix and contrasts
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

desired_comparisons<-treatment
design <- model.matrix(~0+desired_comparisons) 
colnames(design) <- gsub("desired_comparisons", "", colnames(design))

contr.matrix <- makeContrasts(pfc_V_F_vs_pfc_C_F = V - C,
                              levels=colnames(design))



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Remove heteroscedascity from count data
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
v <- voom(x_pfc, design, plot=TRUE)
#write.table(v$E, "~/voom_pfc_F.txt", sep = "\t")

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#  Fit linear models for comparisons of interest
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit, adjust.method="none", p.value = .05)
summary(dt)
# pfc_V_F_vs_pfc_C_F
# Down                  296
# NotSig              14514
# Up                    280

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#  get results
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

results<-topTable(efit, n=Inf, adjust.method = "BH")


