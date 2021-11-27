############################################################################
# ############################################################################
# ######    
# ######    Blood-Brain Correlation Analysis
# ######   
# ######    
# ############################################################################

library(psych)
options(stringsAsFactors = F)
# path to where normalized gene expression data are located
directory<-"~/"
setwd(directory)

# read in normalized gene expression datasets
# these files can be output from the limma voom scripts
  
blood_m<-read.delim(paste0(directory, "voom_bld_M.txt"))
blood_f<-read.delim(paste0(directory, "voom_bld_F.txt"))
amy_m<-read.delim(paste0(directory, "voom_amy_M.txt"))
amy_f<-read.delim(paste0(directory, "voom_amy_F.txt"))
hyp_m<-read.delim(paste0(directory, "voom_hyp_M.txt"))
hyp_f<-read.delim(paste0(directory, "voom_hyp_F.txt"))
pfc_m<-read.delim(paste0(directory, "voom_pfc_M.txt"))
pfc_f<-read.delim(paste0(directory, "voom_pfc_F.txt"))


#~#~#~#~#~#~#~#~#~#~#
#
# Male
# Blood and PFC
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
male_blood_pfc_pop<-intersect(row.names(blood_m), row.names(pfc_m))
blood_m2<-blood_m[row.names(blood_m) %in% male_blood_pfc_pop,]
pfc_m2<-pfc_m[row.names(pfc_m) %in% male_blood_pfc_pop,]
# check they are in the same order
nrow(blood_m2)
#[1] 12908
sum(row.names(blood_m2)==row.names(pfc_m2))
#[1] 12908
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
substr(x = names(blood_m2), start = 1, stop = 2) == substr(x = names(pfc_m2), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_pfc<-list()
for (i in 1:nrow(blood_m2)) {
  
  cor_bld_pfc[[i]]<-corr.test(as.numeric(blood_m2[i,]), as.numeric(pfc_m2[i,]), method = "spearman")
  
}

length(cor_bld_pfc)
#[1] 12908
head(cor_bld_pfc)
names(cor_bld_pfc)<-row.names(blood_m2)
#get only sig correlations
cor_bld_pfc_sig<-Filter(function(y) y[4]<0.05, cor_bld_pfc)
length(cor_bld_pfc_sig)
#[1] 656 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_pfc_sig2<-cor_bld_pfc_sig[list.order(cor_bld_pfc_sig, max(unlist(r)))]
cor_bld_pfc_sig3<-data.frame(ENS=names(cor_bld_pfc_sig2),
                                     corr_spearman = unlist(lapply(cor_bld_pfc_sig2, `[[`, 1)))


#~#~#~#~#~#~#~#~#~#~#
#
# Male
# Blood and AMY
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
male_blood_amy_pop<-intersect(row.names(blood_m), row.names(amy_m))
blood_m2<-blood_m[row.names(blood_m) %in% male_blood_amy_pop,]
amy_m2<-amy_m[row.names(amy_m) %in% male_blood_amy_pop,]
# check they are in the same order
nrow(blood_m2)
#[1] 12935
sum(row.names(blood_m2)==row.names(amy_m2))
#[1] 12935
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
substr(x = names(blood_m2), start = 1, stop = 2) == substr(x = names(amy_m2), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_amy<-list()
for (i in 1:nrow(blood_m2)) {
  
  cor_bld_amy[[i]]<-corr.test(as.numeric(blood_m2[i,]), as.numeric(amy_m2[i,]), method = "spearman")
  
}

length(cor_bld_amy)
#[1] 12935
head(cor_bld_amy)
names(cor_bld_amy)<-row.names(blood_m2)
#get only sig correlations
cor_bld_amy_sig<-Filter(function(y) y[4]<0.05, cor_bld_amy)
length(cor_bld_amy_sig)
#[1] 707 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_amy_sig2<-cor_bld_amy_sig[list.order(cor_bld_amy_sig, max(unlist(r)))]
cor_bld_amy_sig3<-data.frame(ENS=names(cor_bld_amy_sig2),
                             corr_spearman = unlist(lapply(cor_bld_amy_sig2, `[[`, 1)))



#~#~#~#~#~#~#~#~#~#~#
#
# Male
# Blood and hyp
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
male_blood_hyp_pop<-intersect(row.names(blood_m), row.names(hyp_m))
blood_m2<-blood_m[row.names(blood_m) %in% male_blood_hyp_pop,]
hyp_m2<-hyp_m[row.names(hyp_m) %in% male_blood_hyp_pop,]
# check they are in the same order
nrow(blood_m2)
#[1] 12975
sum(row.names(blood_m2)==row.names(hyp_m2))
#[1] 12975
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
substr(x = names(blood_m2), start = 1, stop = 2) == substr(x = names(hyp_m2), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_hyp<-list()
for (i in 1:nrow(blood_m2)) {
  
  cor_bld_hyp[[i]]<-corr.test(as.numeric(blood_m2[i,]), as.numeric(hyp_m2[i,]), method = "spearman")
  
}

length(cor_bld_hyp)
#[1] 12975
head(cor_bld_hyp)
names(cor_bld_hyp)<-row.names(blood_m2)
#get only sig correlations
cor_bld_hyp_sig<-Filter(function(y) y[4]<0.05, cor_bld_hyp)
length(cor_bld_hyp_sig)
#[1] 686 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_hyp_sig2<-cor_bld_hyp_sig[list.order(cor_bld_hyp_sig, max(unlist(r)))]
cor_bld_hyp_sig3<-data.frame(ENS=names(cor_bld_hyp_sig2),
                             corr_spearman = unlist(lapply(cor_bld_hyp_sig2, `[[`, 1)))


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
# Males - Summary of correlated genes
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

bld_amy_corr<-names(cor_bld_amy_sig)
bld_pfc_corr<-names(cor_bld_pfc_sig)
bld_hyp_corr<-names(cor_bld_hyp_sig)

intersect(names(cor_bld_amy_sig),
          names(cor_bld_pfc_sig))

#41

all_corr<-intersect(intersect(names(cor_bld_amy_sig),
                              names(cor_bld_pfc_sig)),
                    names(cor_bld_hyp_sig))

# spearman - 5 transcripts are correlated across blood and each brain area

#[1] "ENSMUSG00000006154.13" "ENSMUSG00000027301.7"  "ENSMUSG00000040565.8"  "ENSMUSG00000090877.3" 
#[5] "ENSMUSG00000091971.3" 

#[1] "Hspa1b" "Hspa1a" "Eps8l1" "Oxt"    "Btaf1" 

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
# Get # positive and negative correlated transcripts
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
bld_amy_corr<-names(cor_bld_amy_sig)
cor_bld_amy_pos_sig<-Filter(function(x) x[1]>0, cor_bld_amy_sig)
length(cor_bld_amy_pos_sig)
#364
cor_bld_amy_neg_sig<-Filter(function(x) x[1]<0, cor_bld_amy_sig)
length(cor_bld_amy_neg_sig)
#343

bld_pfc_corr<-names(cor_bld_pfc_sig)
cor_bld_pfc_pos_sig<-Filter(function(x) x[1]>0, cor_bld_pfc_sig)
length(cor_bld_pfc_pos_sig)
#297
cor_bld_pfc_neg_sig<-Filter(function(x) x[1]<0, cor_bld_pfc_sig)
length(cor_bld_pfc_neg_sig)
#359

bld_hyp_corr<-names(cor_bld_hyp_sig)
cor_bld_hyp_pos_sig<-Filter(function(x) x[1]>0, cor_bld_hyp_sig)
length(cor_bld_hyp_pos_sig)
#283
cor_bld_hyp_neg_sig<-Filter(function(x) x[1]<0, cor_bld_hyp_sig)
length(cor_bld_hyp_neg_sig)
#403


# 1905 unique genes are correlated between blood and at least one brain area
# 113 unique genes are correlated between blood and 2 brain areas
# 5 unique genes are correlated between blood and 3 brain areas

bld_brain_corr[bld_brain_corr$Freq>1,1]
# [1] Btaf1         Eps8l1        Hspa1a        Hspa1b        Oxt           1700017B05Rik 1810055G02Rik
# [8] 4933404O12Rik Abhd10        Adamts10      Agpat5        Akap9         Arl6ip4       Bmpr1b       
# [15] Bok           Bola3         Bptf          Brinp2        C1qc          Ccar1         Ccdc117      
# [22] Chsy1         Cntnap5c      Coasy         Commd5        Cops6         Cox10         Csrnp3       
# [29] Ddx1          Dffa          Dhx36         Dnajb11       Dnajb2        Entpd6        Ess2         
# [36] Fam241b       Foxn2         Gabrg3        Gas1          Gas2l1        Ghr           Gjc3         
# [43] Gm17231       Gm20521       Gm7514        Gsta4         Gtpbp2        Habp4         Hivep1       
# [50] Hypk          Ino80d        Irf2          Jmjd1c        Kdelr1        Kmt5c         Lcp1         
# [57] Llgl1         Magi2         Man2b2        Mbd2          Mgat4c        Mpp3          N4bp2l2      
# [64] Nap1l1        Nbeal1        Nqo2          Ntan1         Nudt7         Nup153        Nupl2        
# [71] Pcnt          Pdia3         Phc3          Ppm1l         Ppp2r2d       Prn           Prpf40b      
# [78] Psenen-ps     Pwp2          Pygl          Rabep1        Rbm34         Rc3h1         Rdh13        
# [85] Rhoc          Rom1          Rps25-ps1     Rragc         Rsrc2         Rwdd2a        Scamp5       
# [92] Scly          Scmh1         Sesn3         Sf3b5         Sgcz          Siva1         Skil         
# [99] Srbd1         Stk3          Tapt1         Tbc1d7        Tdrkh         Tec           Tmem178b     
# [106] Tmem50b       Tnfrsf25      Top1          Trak1         Trp53inp1     Tssc4         Ttll5        
# [113] Ube3b         Vav1          Wdr66         Zdhhc17       Zfp251        Zik1  


#~#~#~#~#~#~#~#~#~#~#
#
# female
# Blood and PFC
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
female_blood_pfc_pop<-intersect(row.names(blood_f), row.names(pfc_f))
blood_f2<-blood_f[row.names(blood_f) %in% female_blood_pfc_pop,]
pfc_f2<-pfc_f[row.names(pfc_f) %in% female_blood_pfc_pop,]
# check they are in the same order
nrow(blood_f2)
#[1] 12226
sum(row.names(blood_f2)==row.names(pfc_f2))
#[1] 12226
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
substr(x = names(blood_f2), start = 1, stop = 2) == substr(x = names(pfc_f2), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_pfc<-list()
for (i in 1:nrow(blood_f2)) {
  
  cor_bld_pfc[[i]]<-corr.test(as.numeric(blood_f2[i,]), as.numeric(pfc_f2[i,]), method = "spearman")
  
}

length(cor_bld_pfc)
#[1] 12226
head(cor_bld_pfc)
names(cor_bld_pfc)<-row.names(blood_f2)
#get only sig correlations
cor_bld_pfc_sig<-Filter(function(y) y[4]<0.05, cor_bld_pfc)
length(cor_bld_pfc_sig)
#[1] 634 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_pfc_sig2<-cor_bld_pfc_sig[list.order(cor_bld_pfc_sig, max(unlist(r)))]
cor_bld_pfc_sig3<-data.frame(ENS=names(cor_bld_pfc_sig2),
                             corr_spearman = unlist(lapply(cor_bld_pfc_sig2, `[[`, 1)))



#~#~#~#~#~#~#~#~#~#~#
#
# female
# Blood and AMY
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
female_blood_amy_pop<-intersect(row.names(blood_f), row.names(amy_f))
blood_f2<-blood_f[row.names(blood_f) %in% female_blood_amy_pop,]
amy_f2<-amy_f[row.names(amy_f) %in% female_blood_amy_pop,]
# check they are in the same order
nrow(blood_f2)
#[1] 12128
sum(row.names(blood_f2)==row.names(amy_f2))
#[1] 12128
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
female_blood_amy_samples<-intersect(substr(names(blood_f), 1, 3), substr(names(amy_f),1,3))
blood_f2_amy<-blood_f2[,substr(names(blood_f2), 1, 3) %in% female_blood_amy_samples]
amy_f2_blood<-amy_f2[,substr(names(amy_f2), 1, 3) %in% female_blood_amy_samples]
substr(x = names(blood_f2_amy), start = 1, stop = 2) == substr(x = names(amy_f2_blood), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_amy<-list()
for (i in 1:nrow(blood_f2)) {
  
  cor_bld_amy[[i]]<-corr.test(as.numeric(blood_f2_amy[i,]), as.numeric(amy_f2_blood[i,]), method = "spearman")
  
}

length(cor_bld_amy)
#[1] 12128
head(cor_bld_amy)
names(cor_bld_amy)<-row.names(blood_f2)
#get only sig correlations
cor_bld_amy_sig<-Filter(function(y) y[4]<0.05, cor_bld_amy)
length(cor_bld_amy_sig)
#[1] 542 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_amy_sig2<-cor_bld_amy_sig[list.order(cor_bld_amy_sig, max(unlist(r)))]
cor_bld_amy_sig3<-data.frame(ENS=names(cor_bld_amy_sig2),
                             corr_spearman = unlist(lapply(cor_bld_amy_sig2, `[[`, 1)))



#~#~#~#~#~#~#~#~#~#~#
#
# female
# Blood and hyp
#
#~#~#~#~#~#~#~#~#~#~#

# only use genes expressed in both blood and pfc
female_blood_hyp_pop<-intersect(row.names(blood_f), row.names(hyp_f))
blood_f2<-blood_f[row.names(blood_f) %in% female_blood_hyp_pop,]
hyp_f2<-hyp_f[row.names(hyp_f) %in% female_blood_hyp_pop,]
# check they are in the same order
nrow(blood_f2)
#[1] 12280
sum(row.names(blood_f2)==row.names(hyp_f2))
#[1] 12280
# ensure that the mouse number order is the same between the two matrices to be corrlated
# so that we indeed get a within subject correlation
female_blood_hyp_samples<-intersect(substr(names(blood_f), 1, 3), substr(names(hyp_f2),1,3))
blood_f2_hyp<-blood_f2[,substr(names(blood_f2), 1, 3) %in% female_blood_hyp_samples]
hyp_f2_blood<-hyp_f2[,substr(names(hyp_f2), 1, 3) %in% female_blood_hyp_samples]
substr(x = names(blood_f2_hyp), start = 1, stop = 2) == substr(x = names(hyp_f2_blood), start = 1, stop = 2)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

cor_bld_hyp<-list()
for (i in 1:nrow(blood_f2)) {
  
  cor_bld_hyp[[i]]<-corr.test(as.numeric(blood_f2_hyp[i,]), as.numeric(hyp_f2_blood[i,]), method = "spearman")
  
}

length(cor_bld_hyp)
#[1] 12280
head(cor_bld_hyp)
names(cor_bld_hyp)<-row.names(blood_f2)
#get only sig correlations
cor_bld_hyp_sig<-Filter(function(y) y[4]<0.05, cor_bld_hyp)
length(cor_bld_hyp_sig)
#[1] 649 spearman
#order sig correlated genes by p value
library(rlist)
cor_bld_hyp_sig2<-cor_bld_hyp_sig[list.order(cor_bld_hyp_sig, max(unlist(r)))]
cor_bld_hyp_sig3<-data.frame(ENS=names(cor_bld_hyp_sig2),
                             corr_spearman = unlist(lapply(cor_bld_hyp_sig2, `[[`, 1)))



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
# females - Summary of correlated genes
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

bld_amy_corr<-names(cor_bld_amy_sig)
bld_pfc_corr<-names(cor_bld_pfc_sig)
bld_hyp_corr<-names(cor_bld_hyp_sig)

intersect(names(cor_bld_amy_sig),
          names(cor_bld_pfc_sig))



all_corr<-intersect(intersect(names(cor_bld_amy_sig),
                              names(cor_bld_pfc_sig)),
                    names(cor_bld_hyp_sig))
length(all_corr)
# spearman - 6 transcripts are correlated across blood and each brain area

#[1] "ENSMUSG00000021250.13" "ENSMUSG00000021270.13" "ENSMUSG00000021548.10" "ENSMUSG00000028782.14"
#[5] "ENSMUSG00000035042.2"  "ENSMUSG00000045314.5" 

DEGs[DEGs$X %in% all_corr, 4]
#[1] "Ccl5"     "Hsp90aa1" "Sowahb"   "Fos"      "Adgrb2"   "Ccnh"    

# 1690 unique genes are correlated between blood and at least one brain area
# 99 unique genes are correlated between blood and 2 brain areas
# 6 unique gene is correlated between blood and 3 brain areas

bld_brain_corr[bld_brain_corr$Freq>1,1]
# [1] Adgrb2        Ccl5          Ccnh          Fos           Hsp90aa1      Sowahb        1810026B05Rik
# [8] 2810002D19Rik 6720489N17Rik Adck1         Ap1b1         Arpc2         Atp13a1       B3gntl1      
# [15] Bcl2l13       Bcl9l         Bin2          Btaf1         Cep83os       Chd7          Clns1a       
# [22] Cnpy4         Col4a3bp      Dipk1b        Doc2b         Dubr          Dusp7         Eef1akmt1    
# [29] Eif3m         Eml3          Erbb4         Fam189b       Fth1          Gabra5        Gars         
# [36] Gfra2         Gm13092       Gm15478       Gm32031       Gm42984       Gm47283       Gm5526       
# [43] Gm8702        Gon4l         Hmox2         Hnrnpf        Hnrnph1       Hs3st4        Hsd17b4      
# [50] Ifnar2        Inpp1         Iqsec1        Kcnk1         Kctd3         Lsg1          Map3k14      
# [57] Map3k6        Mettl3        Mrfap1        Mrpl11        Mrpl39        Mrpl45        Myof         
# [64] Nfkbib        Nudt7         Pabpc1        Pias3         Plin2         Plvap         Ppp4r2       
# [71] Prpsap1       Psmb2         Ptch1         Pum1          Rab5b         Rabgap1       Rabgap1l     
# [78] Rad54b        Rbbp4         Rbm27         Resp18        Rny1          Sae1          Sec24d       
# [85] Sema5b        Sh2b2         Slc12a5       Slc25a35      Slc7a6        Snx8          Srfbp1       
# [92] Tmem123       Tmem131       Tmem138       Tmx3          Tom1          Tspan31       Ttr          
# [99] Usp54         Wdr6          Wsb1          Xbp1          Zc3h14        Zfp369        Zfp444  

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
# Get # positive and negative correlated transcripts
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
bld_amy_corr<-names(cor_bld_amy_sig)
cor_bld_amy_pos_sig<-Filter(function(x) x[1]>0, cor_bld_amy_sig)
length(cor_bld_amy_pos_sig)
#229
cor_bld_amy_neg_sig<-Filter(function(x) x[1]<0, cor_bld_amy_sig)
length(cor_bld_amy_neg_sig)
#313

bld_pfc_corr<-names(cor_bld_pfc_sig)
cor_bld_pfc_pos_sig<-Filter(function(x) x[1]>0, cor_bld_pfc_sig)
length(cor_bld_pfc_pos_sig)
#367
cor_bld_pfc_neg_sig<-Filter(function(x) x[1]<0, cor_bld_pfc_sig)
length(cor_bld_pfc_neg_sig)
#267

bld_hyp_corr<-names(cor_bld_hyp_sig)
cor_bld_hyp_pos_sig<-Filter(function(x) x[1]>0, cor_bld_hyp_sig)
length(cor_bld_hyp_pos_sig)
#305
cor_bld_hyp_neg_sig<-Filter(function(x) x[1]<0, cor_bld_hyp_sig)
length(cor_bld_hyp_neg_sig)
#344
