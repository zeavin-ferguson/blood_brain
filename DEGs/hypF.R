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

hyp_files<-list.files(directory)[grep(pattern = "F.*HYP|F.*hyp", x = list.files(directory))]

x_hyp <- readDGE(hyp_files, columns= c(1,2)) 
class(x_hyp)
dim(x_hyp)
#[1] 55481    18

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize sample information
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Assign sample-level information related to the experimental design to the columns of the counts matrix [treatment group (CIE vs Air), sex (M and F), tissue (blood, hypgdala, prefrontal cortex, hypothalamus)] 
samplenames <- substring(colnames(x_hyp), 1, 21)
colnames(x_hyp) <- samplenames
treatment <- as.factor(substr(samplenames, 17, 17))

x_hyp$samples$treatment <- treatment

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#       Organize gene annotations
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
geneid <- rownames(x_hyp) 
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
genes<-genes[match(x = gsub("(.*)(\\..*)","\\1",rownames(x_hyp)), table = genes$ensembl_gene_id),]
x_hyp$genes<-genes

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      remove ERCC counts from the count files
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

ercc.exprs <- grep(pattern = "GERCC", x = rownames(x_hyp))
# remove ERCCs
x_hyp <-x_hyp[-ercc.exprs,, keep.lib.sizes=FALSE]
dim(x_hyp)
## [1] 55389    18


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#      Remove genes that are lowly expressed
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

keep.exprs <- rowSums(cpm(x_hyp)>1)>=8
x_hyp <- x_hyp[keep.exprs,, keep.lib.sizes=FALSE]
dim(x_hyp)
#[1] 15469   18

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Normalize gene expression distributions
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

x_hyp <- calcNormFactors(x_hyp, method = "TMM")

x_hyp$samples$norm.factors
# [1] 0.9661840 0.9601353 0.9586465 0.9942822 1.0163999 0.9718565 0.9684410 0.9870274 0.9985593
# [10] 1.0238279 1.0467824 1.0031912 1.0185839 0.9994081 1.0335779 1.0079814 1.0239594 1.0273333

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

contr.matrix <- makeContrasts(hyp_V_F_vs_hyp_C_F = V - C,
                              levels=colnames(design))



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#   Remove heteroscedascity from count data
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
v <- voom(x_hyp, design, plot=TRUE)
#write.table(v$E, "~/voom_hyp_F.txt", sep = "\t")

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
# hyp_V_F_vs_hyp_C_F
# Down                  386
# NotSig              14642
# Up                    441

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
#  get results
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

results<-topTable(efit, n=Inf, adjust.method = "BH")


