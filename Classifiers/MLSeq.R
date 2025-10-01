#MLSeq
#https://rdrr.io/bioc/MLSeq/f/inst/doc/MLSeq.pdf

# This script runs MLSeq on the blood count data to predict CIE vs Control samples

# MLSeq was created specifically for applying machine learning classification to RNASeq data
# it uses the caret package for model training and tuning
# one major benefit of the MLSeq package is that it takes in raw count files as the input and includes the normalization and transformation steps in the model training
# then it applies the same parameters learned during training to the test samples
# this removes a source of data leakage which could occur if information from the test samples is used for normalization and transformation
# this would not be a problem for log transformations but is a problem is gene variance across all samples is used for shrinkage for example

# trained(fit) shows cross validation results of tuning parameters
# plot(fit) plots the CV results
# varImp(trained(fit)) gives variable importance.

# load libraries
library(MLSeq)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(pamr)
library(EnhancedVolcano)
library(pROC)
library(R.utils)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Get raw count files from Gene Expression Omnibus (GEO)
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# go to this GEO link and download the raw count files (GSE176122_RAW.tar) listed under Supplementary Files at the bottom of the page
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176122

# untar the file into the count_directory
count_directory<-"PATH/TO/RAW/COUNT/FILES/"

count_directory<-"~/Dropbox/blood_brain_project/TagSeq_Analysis/trimmed_data/Bowtie2/Htseq_Bowtie2/"
count_directory<-"~/Downloads/GSE176122_RAW (1)"

# List all .gz files in the directory
files <- list.files(path = count_directory, pattern = "\\.gz$", full.names = TRUE)

# unzip .gz files in the count_directory if not already done
for (f in files) {
  gunzip(f, overwrite = F)
}

# get raw count files in the directory for blood samples
# do this for only females, only males, and both sexes so we can build models for each

#female blood samples
sampleFiles1<-list.files(count_directory)[grep(pattern = "F-.*bld|F-.*BLD", x = list.files(count_directory))]
#male blood samples
sampleFiles2<-list.files(count_directory)[grep(pattern = "M-.*bld|M-.*BLD", x = list.files(count_directory))]
#all blood samples
sampleFiles3<-list.files(count_directory)[grep(pattern = "bld|BLD", x = list.files(count_directory))]
sampleFiles<-list(bloodF=sampleFiles1,
                  bloodM=sampleFiles2,
                  blood=sampleFiles3)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Create DESeqDataSets from raw count files
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# MLSeq requires the input to be a DESeqDataSet object

dds<-list()
for (i in 1:3) {
  temp<-toupper(gsub("(GSM.*_)(.*_S.*)(_R1.*)", "\\2", sampleFiles[[i]]))
  sampleCondition <- factor(substr(x = temp, start = 6, stop = 6))
  sampleTable <- data.frame(sampleName = substr(x = temp, start = 1, stop = 10),
                            fileName = sampleFiles[[i]],
                            condition = sampleCondition)
  dds[[i]] <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = count_directory,
                                         design= ~ condition)
}

names(dds) <- names(sampleFiles)

lapply(dds, function(x) head(counts(x)))
lapply(dds, function(x) dim(counts(x)))
# $bloodF
# [1] 55477    19
# 
# $bloodM
# [1] 55477    18
# 
# $blood
# [1] 55477    37

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Filter datasets to remove low count genes and ERCC spike-ins
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# remove ERCC
dds<-lapply(dds, function(x) x[-grep(x=row.names(counts(x)), pattern = "GERCC"),])
keep <- rowSums(counts(dds[[i]])) >= ceiling(ncol(dds[[i]])/2)
dds <- lapply(dds, function(x) x[rowSums(counts(x)) >= ceiling(ncol(x)/2),])
lapply(dds, dim)
# $bloodF
# [1] 15823    19
# 
# $bloodM
# [1] 17071    18
# 
# $blood
# [1] 16484    37

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Build models
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

fit.rf<-list()
fit.pls<-list()
fit.glmnet<-list()

for (i in 1:length(dds)) {

  data <- counts(dds[[i]])
  classtr <- DataFrame(condition = as.factor(dds[[i]]$condition))
  # I will use repeated CV to estimate accuracy
  # Minimum count is set to 1 in order to prevent 0 division problem within classification models.
  data.train <- as.matrix(data + 1)

  data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                        design = formula(~condition))
   # Define control lists.

  set.seed(1)
  ctrl.continuous <- trainControl(method = "repeatedcv", 
                                  number = 5 , 
                                  repeats = 10, 
                                  classProbs = TRUE,  # must include classProb=T to make ROC curves later
                                  savePredictions = "final")
  
  # Continuous classifiers

  set.seed(1)
  fit.rf[[i]] <- classify(data = data.trainS4, method = "rf",
                          preProcessing = "deseq-vst", ref = "C", tuneLength = 9,
                          control = ctrl.continuous)

  set.seed(1)
  fit.pls[[i]] <- classify(data = data.trainS4, method = "pls",
                           preProcessing = "deseq-vst", ref = "C", tuneLength = 9,
                           control = ctrl.continuous)

  set.seed(1)
  fit.glmnet[[i]] <- classify(data = data.trainS4, method = "glmnet",
                              preProcessing = "deseq-vst", ref = "C", tuneLength = 9,
                              control = ctrl.continuous)

}

names(fit.rf)<-names(dds)
names(fit.pls)<-names(dds)
names(fit.glmnet)<-names(dds)

# create a directory to save the models
MLSeq_results <- "PATH/TO/RESULT/FILES/"

# save dds and models
save(dds, fit.rf, fit.pls, fit.glmnet, file = paste0(MLSeq_results, "dds_and_models.save"))

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Plot ROC curves for each model
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# LOAD MODELS AND DDS
load(file = paste0(MLSeq_results, "dds_and_models.save"))

Males <- list(RF=fit.rf$bloodM,
              LR=fit.glmnet$bloodM,
              PLSDA=fit.pls$bloodM)

Females <- list(RF=fit.rf$bloodF,
                LR=fit.glmnet$bloodF,
                PLSDA=fit.pls$bloodF)

malesANDfemales <- list(RF=fit.rf$blood,
                        LR=fit.glmnet$blood,
                        PLSDA=fit.pls$blood)

ROC_colors <- c("blue", "red", "black")
AUC_print_local <- c(0.4, 0.35, 0.3)*100

pdf(file = paste0(MLSeq_results, "ROC_blood.pdf"),
    width = 10, height = 5)

par(pty="s") 

roc(malesANDfemales[[1]]@modelInfo@trainedModel$pred$obs, 
    malesANDfemales[[1]]@modelInfo@trainedModel$pred$V, 
    plot=TRUE, 
    print.auc=T, 
    col=ROC_colors[[1]], 
    print.auc.y=AUC_print_local[[1]],
    print.auc.cex=1,
    lwd = 3, 
    legacy.axes=TRUE, 
    percent = TRUE,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    main="Peripheral Blood")

for (i in 2:length(malesANDfemales)) {
  
  roc(malesANDfemales[[i]]@modelInfo@trainedModel$pred$obs, 
      malesANDfemales[[i]]@modelInfo@trainedModel$pred$V, 
      plot=TRUE, 
      print.auc=T, 
      col=ROC_colors[[i]], 
      print.auc.y=AUC_print_local[[i]],
      print.auc.cex=1,
      lwd = 3, 
      legacy.axes=TRUE,
      percent = TRUE,
      add = TRUE)
  
}

legend("bottomright",legend=names(malesANDfemales),col=ROC_colors[1:length(malesANDfemales)],lwd=3, cex = 0.75)

dev.off()



ROC_colors <- c("blue", "red", "black")
AUC_print_local <- c(0.4, 0.35, 0.3)*100

pdf(file = paste0(MLSeq_results, "ROC_blood_M.pdf"),
    width = 10, height = 5)

par(pty="s") 

roc(Males[[1]]@modelInfo@trainedModel$pred$obs, 
    Males[[1]]@modelInfo@trainedModel$pred$V, 
    plot=TRUE, 
    print.auc=T, 
    col=ROC_colors[[1]], 
    print.auc.y=AUC_print_local[[1]],
    print.auc.cex=1,
    lwd = 3, 
    legacy.axes=TRUE, 
    percent = TRUE,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    main="Male Peripheral Blood")

for (i in 2:length(Males)) {
  
  roc(Males[[i]]@modelInfo@trainedModel$pred$obs, 
      Males[[i]]@modelInfo@trainedModel$pred$V, 
      plot=TRUE, 
      print.auc=T, 
      col=ROC_colors[[i]], 
      print.auc.y=AUC_print_local[[i]],
      print.auc.cex=1,
      lwd = 3, 
      legacy.axes=TRUE,
      percent = TRUE,
      add = TRUE)
  
}

legend("bottomright",legend=names(Males),col=ROC_colors[1:length(Males)],lwd=3, cex = 0.75)

dev.off()




ROC_colors <- c("blue", "red", "black")
AUC_print_local <- c(0.4, 0.35, 0.3)*100

pdf(file = paste0(MLSeq_results, "ROC_blood_F.pdf"),
    width = 10, height = 5)
par(pty="s") 

roc(Females[[1]]@modelInfo@trainedModel$pred$obs, 
    Females[[1]]@modelInfo@trainedModel$pred$V, 
    plot=TRUE, 
    print.auc=T, 
    col=ROC_colors[[1]], 
    print.auc.y=AUC_print_local[[1]],
    print.auc.cex=1,
    lwd = 3, 
    legacy.axes=TRUE, 
    percent = TRUE,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    main="Female Peripheral Blood")

for (i in 2:length(Females)) {
  
  roc(Females[[i]]@modelInfo@trainedModel$pred$obs, 
      Females[[i]]@modelInfo@trainedModel$pred$V, 
      plot=TRUE, 
      print.auc=T, 
      col=ROC_colors[[i]], 
      print.auc.y=AUC_print_local[[i]],
      print.auc.cex=1,
      lwd = 3, 
      legacy.axes=TRUE,
      percent = TRUE,
      add = TRUE)
  
}

legend("bottomright",legend=names(Females),col=ROC_colors[1:length(Females)],lwd=3, cex = 0.75)

dev.off()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# 
# Extract variable importance from models
# 
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Extract variable importance
# overall_model_features <- varImp(trained(fit.glmnet$blood))$importance
# overall_model_features <- varImp(trained(fit.glmnet$bloodF))$importance
overall_model_features <- varImp(trained(fit.glmnet$bloodM))$importance

# Order by importance, decreasing
ordered_features <- overall_model_features[order(overall_model_features$Overall, 
                                                 decreasing = TRUE), , drop = FALSE]

# Keep rownames (these are your Ensembl IDs)
head(ordered_features)   # top features

nonzero_features <- ordered_features[ordered_features$Overall > 0, , drop = FALSE]
dim(nonzero_features)   # should show 39 x 1
rownames(nonzero_features)   # Ensembl IDs of those features

nonzero_features_df <- data.frame(
  ensembl_id = rownames(nonzero_features),
  importance = nonzero_features$Overall,
  row.names = NULL
)
head(nonzero_features_df)

library(AnnotationDbi)
library(org.Mm.eg.db)

# Strip version numbers from Ensembl IDs
nonzero_features_df$ensembl_clean <- sub("\\..*", "", nonzero_features_df$ensembl_id)

# Map Ensembl IDs to gene symbols using org.Mm.eg.db
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = nonzero_features_df$ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"  # if multiple mappings, take the first
)

# Add gene symbols back into your data frame
nonzero_features_df$gene_symbol <- gene_symbols[nonzero_features_df$ensembl_clean]

# Order by importance
nonzero_features_df <- nonzero_features_df[order(nonzero_features_df$importance, 
                                                 decreasing = TRUE), ]

nonzero_features_df

