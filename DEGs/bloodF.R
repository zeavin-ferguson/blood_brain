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

bld_files<-list.files(directory)[grep(pattern = "F.*bld|F.*bld", x = list.files(directory))]

x_bld <- readDGE(bld_files, columns= c(1,2)) 
class(x_bld)
dim(x_bld)
#[1] 55481    19

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize sample information
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Assign sample-level information related to the experimental design to the columns of the counts matrix [treatment group (CIE vs Air), sex (M and F), tissue (blood, bldgdala, prefrontal cortex, hypothalamus)] 
samplenames <- substring(colnames(x_bld), 1, 21)
colnames(x_bld) <- samplenames
treatment <- as.factor(substr(samplenames, 17, 17))

x_bld$samples$treatment <- treatment

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize gene annotations
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
geneid <- rownames(x_bld) 
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
genes<-genes[match(x = gsub("(.*)(\\..*)","\\1",rownames(x_bld)), table = genes$ensembl_gene_id),]
x_bld$genes<-genes

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      remove ERCC counts from the count files
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

ercc.exprs <- grep(pattern = "GERCC", x = rownames(x_bld))
# remove ERCCs
x_bld <-x_bld[-ercc.exprs,, keep.lib.sizes=FALSE]
dim(x_bld)
## [1] 55389    19


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      Remove genes that are lowly expressed
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

keep.exprs <- rowSums(cpm(x_bld)>1)>=8
x_bld <- x_bld[keep.exprs,, keep.lib.sizes=FALSE]
dim(x_bld)
#[1] 13380   19

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Normalize gene expression distributions
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

x_bld <- calcNormFactors(x_bld, method = "TMM")

x_bld$samples$norm.factors
#[1] 1.0582444 1.5707113 1.0628610 0.7969108 1.1536923 1.7276392 1.5275455 1.0037808
#[9] 0.9729211 1.1492486 0.7857360 0.7920940 0.7242358 0.7067661 0.9313472 0.8247515
#[17] 0.8819361 0.8022910 1.2004690

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

contr.matrix <- makeContrasts(bld_V_F_vs_bld_C_F = V - C,
                              levels=colnames(design))



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Remove heteroscedascity from count data
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
v <- voom(x_bld, design, plot=TRUE)
#write.table(v$E, "~/voom_bld_F.txt", sep = "\t")

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
# bld_V_F_vs_bld_C_F
# Down                  709
# NotSig              11814
# Up                    857

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#  get results
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

results<-topTable(efit, n=Inf, adjust.method = "BH")




