
#set options
options(stringsAsFactors = F)

#load libraries
library(ggplot2)
library(psych)
library(cocor)
library(AnnotationDbi)
library(org.Mm.eg.db)

# path to where normalized gene expression data are located
directory<-"~/"
setwd(directory)

# read in normalized expression data

exprAll<-list(BLD_m=read.delim(paste0(directory, "voom_bld_M.txt")),
              BLD_f=read.delim(paste0(directory, "voom_bld_F.txt")),
              AMY_m=read.delim(paste0(directory, "voom_amy_M.txt")),
              AMY_f=read.delim(paste0(directory, "voom_amy_F.txt")),
              HYP_m=read.delim(paste0(directory, "voom_hyp_M.txt")),
              HYP_f=read.delim(paste0(directory, "voom_hyp_F.txt")),
              PFC_m=read.delim(paste0(directory, "voom_pfc_M.txt")),
              PFC_f=read.delim(paste0(directory, "voom_pfc_F.txt")))

lapply(exprAll, dim)

# get average gene expression values and the standard errors for 
# (1) all subjects, (2) only ethanol vapor animals, and (3) just air vapor animals

std <- function(x) sd(x)/sqrt(length(x))

for (i in 1:length(exprAll)) {
  exprAll[[i]]$avgExpr<-rowMeans(exprAll[[i]][,grepl(pattern = "_V_|_C_", x = names(exprAll[[i]]))])
  exprAll[[i]]$avgExprAlcDependent<-rowMeans(exprAll[[i]][,grepl(pattern = "_V_", x = names(exprAll[[i]]))])
  exprAll[[i]]$avgExprControl<-rowMeans(exprAll[[i]][,grepl(pattern = "_C_", x = names(exprAll[[i]]))])
  
  exprAll[[i]]$avgExprSE<-apply(exprAll[[i]][,grepl(pattern = "_V_|_C_", x = names(exprAll[[i]]))], MARGIN = 1, std)
  exprAll[[i]]$avgExprAlcDependentSE<-apply(exprAll[[i]][,grepl(pattern = "_V_", x = names(exprAll[[i]]))], MARGIN = 1, std)
  exprAll[[i]]$avgExprControlSE<-apply(exprAll[[i]][,grepl(pattern = "_C_", x = names(exprAll[[i]]))], MARGIN = 1, std)
}


# add gene info to expression datasets (this is so that we can subset the dataframes by cell type specific genes or interferon response genes because those gene lists are gene symbols and not mouse emsembl ids)

for (i in 1:length(exprAll)) {
  exprAll[[i]]$Symbol<-toupper(mapIds(x = org.Mm.eg.db, keys = gsub(pattern = "(.*)(\\..*)", replacement = "\\1", x = row.names(exprAll[[i]])), column = "SYMBOL", keytype = "ENSEMBL"))
}

# extract the last 7 columns which are the summary info and gene symbols just added
exprAllSummary<-lapply(exprAll, function(x) x[,(ncol(x)-6):ncol(x)])
str(exprAllSummary)
lapply(exprAllSummary, dim)

# append the list element name to each column in the list element data frame

for (i in 1:length(exprAllSummary)) {
  names(exprAllSummary[[i]])<-paste0(names(exprAllSummary[[i]]), ".", names(exprAllSummary)[[i]])
}

exprCorsMales<-list(AMY=merge(exprAllSummary$BLD_m, exprAllSummary$AMY_m, by="row.names"),
                    PFC=merge(exprAllSummary$BLD_m, exprAllSummary$PFC_m, by="row.names"),
                    HYP=merge(exprAllSummary$BLD_m, exprAllSummary$HYP_m, by="row.names"))

exprCorsFemales<-list(AMY=merge(exprAllSummary$BLD_f, exprAllSummary$AMY_f, by="row.names"),
                      PFC=merge(exprAllSummary$BLD_f, exprAllSummary$PFC_f, by="row.names"),
                      HYP=merge(exprAllSummary$BLD_f, exprAllSummary$HYP_f, by="row.names"))

# remove the first 4 rows which are __ambiguous, __no_feature, __not_aligned, and __too_low_aQual
exprCorsMales<-lapply(exprCorsMales, function(x) x[-c(1:4),])
exprCorsFemales<-lapply(exprCorsFemales, function(x) x[-c(1:4),])

# add cell type response info
library(BRETIGEA)
str(markers_df_brain)
unique(markers_df_brain$cell)

library(imsig)
sig
unique(sig$cell)
sig[sig$cell=="Interferon",]

for (i in 1:length(exprCorsMales)) {
  exprCorsMales[[i]]<-merge(exprCorsMales[[i]], sig, by.x="Symbol.BLD_m", by.y="gene", all.x=T)
  exprCorsMales[[i]]<-merge(exprCorsMales[[i]], markers_df_brain, by.x="Symbol.BLD_m", by.y="markers", all.x=T)
}

for (i in 1:length(exprCorsFemales)) {
  exprCorsFemales[[i]]<-merge(exprCorsFemales[[i]], sig, by.x="Symbol.BLD_f", by.y="gene", all.x=T)
  exprCorsFemales[[i]]<-merge(exprCorsFemales[[i]], markers_df_brain, by.x="Symbol.BLD_f", by.y="markers", all.x=T)
}


# compute correlations

corsMales<-list()
corsFemales<-list()
corsMalesControl<-list()
corsMalesVapor<-list()
corsFemalesControl<-list()
corsFemalesVapor<-list()

names(exprCorsMales[[1]])
# column 3: "avgExpr.BLD_m"
# column 4: "avgExprAlcDependent.BLD_m"
# column 5: "avgExprControl.BLD_m" 
# column 9: "avgExpr.HYP_m"
# column 10: "avgExprAlcDependent.HYP_m"
# column 11: "avgExprControl.HYP_m"


for (i in 1:length(exprCorsMales)) {
  corsMales[[i]]<-corr.test(exprCorsMales[[i]][,3], exprCorsMales[[i]][,9], method = "spearman")
  corsMalesControl[[i]]<-corr.test(exprCorsMales[[i]][,5], exprCorsMales[[i]][,11], method = "spearman")
  corsMalesVapor[[i]]<-corr.test(exprCorsMales[[i]][,4], exprCorsMales[[i]][,10], method = "spearman")
  corsFemales[[i]]<-corr.test(exprCorsFemales[[i]][,3], exprCorsFemales[[i]][,9], method = "spearman")
  corsFemalesControl[[i]]<-corr.test(exprCorsFemales[[i]][,5], exprCorsFemales[[i]][,11], method = "spearman")
  corsFemalesVapor[[i]]<-corr.test(exprCorsFemales[[i]][,4], exprCorsFemales[[i]][,10], method = "spearman")
}

names(corsMales)<-names(exprCorsMales)
names(corsMalesControl)<-names(exprCorsMales)
names(corsMalesVapor)<-names(exprCorsMales)
names(corsFemales)<-names(exprCorsFemales)
names(corsFemalesControl)<-names(exprCorsFemales)
names(corsFemalesVapor)<-names(exprCorsFemales)

pdf("~/correlation_between_subjects_for_figs_male.pdf")
for (i in 1:length(exprCorsMales)) {
  p<-ggplot(data = exprCorsMales[[i]], mapping = aes(x = exprCorsMales[[i]][,3], y=exprCorsMales[[i]][,9])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    #geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nMale Mice") +
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExpr\\.)(.*)(_m)", replacement = "\\2", x = names(exprCorsMales[[i]])[9]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, size=8, label=paste("rho =",signif(corsMales[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()

pdf("~/correlation_between_subjects_for_figs_female.pdf")
for (i in 1:length(exprCorsFemales)) {
  p<-ggplot(data = exprCorsFemales[[i]], mapping = aes(x = exprCorsFemales[[i]][,3], y=exprCorsFemales[[i]][,9])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    #geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nFemale Mice") +
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExpr\\.)(.*)(_f)", replacement = "\\2", x = names(exprCorsFemales[[i]])[9]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, size=8, label=paste("rho =",signif(corsFemales[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()


cor_comparison_results<-list()
for (i in 1:length(exprCorsFemales)) {
  cor_comparison_results[[i]]<-cocor.indep.groups(r1.jk=corsFemales[[i]]$r, r2.hm=corsMales[[i]]$r, n1=corsFemales[[i]]$n, n2=corsMales[[i]]$n, alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)
}
names(cor_comparison_results)<-names(exprCorsFemales)
# $AMY
# 
# Results of a comparison of two correlations based on independent groups
# 
# Comparison between r1.jk = 0.5051 and r2.hm = 0.6729
# Difference: r1.jk - r2.hm = -0.1678
# Group sizes: n1 = 12242, n2 = 13062
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
# 
# fisher1925: Fisher's z (1925)
#   z = -20.6563, p-value = 0.0000
#   Null hypothesis rejected
# 
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.1840 -0.1516
# Null hypothesis rejected (Interval does not include 0)
# 
# 
# $PFC
# 
# Results of a comparison of two correlations based on independent groups
# 
# Comparison between r1.jk = 0.4996 and r2.hm = 0.6668
# Difference: r1.jk - r2.hm = -0.1672
# Group sizes: n1 = 12340, n2 = 13036
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
# 
# fisher1925: Fisher's z (1925)
#   z = -20.3955, p-value = 0.0000
#   Null hypothesis rejected
# 
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.1836 -0.1509
# Null hypothesis rejected (Interval does not include 0)
# 
# 
# $HYP
# 
# Results of a comparison of two correlations based on independent groups
# 
# Comparison between r1.jk = 0.4959 and r2.hm = 0.6667
# Difference: r1.jk - r2.hm = -0.1708
# Group sizes: n1 = 12395, n2 = 13102
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
# 
# fisher1925: Fisher's z (1925)
#   z = -20.8275, p-value = 0.0000
#   Null hypothesis rejected
# 
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.1872 -0.1545
# Null hypothesis rejected (Interval does not include 0)












for (i in 1:length(exprCorsMales)) {
  corsMalesControl[[i]]<-corr.test(exprCorsMales[[i]][,5], exprCorsMales[[i]][,11], method = "spearman")
  corsMalesVapor[[i]]<-corr.test(exprCorsMales[[i]][,4], exprCorsMales[[i]][,10], method = "spearman")
  corsFemalesControl[[i]]<-corr.test(exprCorsFemales[[i]][,5], exprCorsFemales[[i]][,11], method = "spearman")
  corsFemalesVapor[[i]]<-corr.test(exprCorsFemales[[i]][,4], exprCorsFemales[[i]][,10], method = "spearman")
}

names(corsMalesControl)<-names(exprCorsMales)
names(corsMalesVapor)<-names(exprCorsMales)
names(corsFemalesControl)<-names(exprCorsFemales)
names(corsFemalesVapor)<-names(exprCorsFemales)

pdf("~/correlation_between_subjects_for_figs_male_control.pdf")
for (i in 1:length(exprCorsMales)) {
  p<-ggplot(data = exprCorsMales[[i]], mapping = aes(x = exprCorsMales[[i]][,5], y=exprCorsMales[[i]][,11])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nNon-Dependent Male Mice") +
    #xlab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsMales[[i]])[5])) +
    #ylab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsMales[[i]])[11])) 
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExprControl\\.)(.*)(_m)", replacement = "\\2", x = names(exprCorsMales[[i]])[11]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, label=paste("rho =",signif(corsMalesControl[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()



pdf("~/correlation_between_subjects_for_figs_male_alcohol_vapor.pdf")
for (i in 1:length(exprCorsMales)) {
  p<-ggplot(data = exprCorsMales[[i]], mapping = aes(x = exprCorsMales[[i]][,4], y=exprCorsMales[[i]][,10])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nAlcohol-Dependent Male Mice") +
    #xlab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsMales[[i]])[5])) +
    #ylab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsMales[[i]])[11])) 
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExprAlcDependent\\.)(.*)(_m)", replacement = "\\2", x = names(exprCorsMales[[i]])[10]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, label=paste("rho =",signif(corsMalesVapor[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()





pdf("~/correlation_between_subjects_for_figs_female_control.pdf")
for (i in 1:length(exprCorsFemales)) {
  p<-ggplot(data = exprCorsFemales[[i]], mapping = aes(x = exprCorsFemales[[i]][,5], y=exprCorsFemales[[i]][,11])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nNon-Dependent Female Mice") +
    #xlab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsFemales[[i]])[5])) +
    #ylab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsFemales[[i]])[11])) 
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExprControl\\.)(.*)(_f)", replacement = "\\2", x = names(exprCorsFemales[[i]])[11]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, label=paste("rho =",signif(corsFemalesControl[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()



pdf("~/correlation_between_subjects_for_figs_female_alcohol_vapor.pdf")
for (i in 1:length(exprCorsFemales)) {
  p<-ggplot(data = exprCorsFemales[[i]], mapping = aes(x = exprCorsFemales[[i]][,4], y=exprCorsFemales[[i]][,10])) +
    geom_point(size = .2)+
    xlim(c(0,13))+
    ylim(c(0,13))+
    geom_smooth(method="lm")+
    ggtitle("Blood-Brain Correlations\nAlcohol-Dependent Female Mice") +
    #xlab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsFemales[[i]])[5])) +
    #ylab(gsub(pattern = "(avgExpr)(.*)", replacement = "\\2", x = names(exprCorsFemales[[i]])[11])) 
    xlab("Blood Normalized Expression") +
    ylab(paste(gsub(pattern = "(avgExprAlcDependent\\.)(.*)(_f)", replacement = "\\2", x = names(exprCorsFemales[[i]])[10]), "Normalized Expression")) +
    geom_text(x=2.5, y=11, label=paste("rho =",signif(corsFemalesVapor[[i]]$r, 2)))+
    theme(plot.title = element_text(color="black", size=14, face="bold.italic"),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"))
  print(p)
}
dev.off()

