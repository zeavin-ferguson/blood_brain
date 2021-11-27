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

amy_files<-list.files(directory)[grep(pattern = "M-.*amy|M-.*AMY", x = list.files(directory))]

x_amy <- readDGE(amy_files, columns= c(1,2)) 
class(x_amy)
dim(x_amy)
#[1] 55481    18

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize sample information
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Assign sample-level information related to the experimental design to the columns of the counts matrix [treatment group (CIE vs Air), sex (M and F), tissue (blood, amygdala, prefrontal cortex, hypothalamus)] 
samplenames <- substring(colnames(x_amy), 1, 21)
colnames(x_amy) <- samplenames
treatment <- as.factor(substr(samplenames, 17, 17))

x_amy$samples$treatment <- treatment

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize gene annotations
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
geneid <- rownames(x_amy) 
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
genes<-genes[match(x = gsub("(.*)(\\..*)","\\1",rownames(x_amy)), table = genes$ensembl_gene_id),]
x_amy$genes<-genes

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      remove ERCC counts from the count files
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

ercc.exprs <- grep(pattern = "GERCC", x = rownames(x_amy))
# remove ERCCs
x_amy <-x_amy[-ercc.exprs,, keep.lib.sizes=FALSE]
dim(x_amy)
## [1] 55389    18


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      Remove genes that are lowly expressed
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

keep.exprs <- rowSums(cpm(x_amy)>1)>=8
x_amy <- x_amy[keep.exprs,, keep.lib.sizes=FALSE]
dim(x_amy)
#[1] 14961   18

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Normalize gene expression distributions
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

x_amy <- calcNormFactors(x_amy, method = "TMM")

x_amy$samples$norm.factors
#[1] 0.9693862 0.9789791 1.0068001 0.9856157 0.9753882 0.9978552 1.0371960 1.0304524
#[9] 1.0003078 0.9969655 1.0084451 1.0240327 1.0064559 0.9905849 1.0006538 1.0152051
#[17] 0.9876433 0.9909256

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

contr.matrix <- makeContrasts(amy_V_M_vs_amy_C_M = V - C,
                              levels=colnames(design))



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Remove heteroscedascity from count data
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
v <- voom(x_amy, design, plot=TRUE)
#write.table(v$E, "~/voom_amy_M.txt", sep = "\t")

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
# amy_V_M_vs_amy_C_M
# Down                  626
# NotSig              13628
# Up                    707

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#  get results
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

results<-topTable(efit, n=Inf, adjust.method = "BH")
