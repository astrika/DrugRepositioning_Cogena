# Astrid M Manuel
# Runing Cogena for drug repositioning strategies
# Last Update 04/05/2021

library(cogena)
setwd("C:/Users/amanuel1/Dev/ZhaoLab/Cogena")
RNAseq <- read.table("C:/users/amanuel1/Dev/ZhaoLab/MS/GSE111972_norm_data.txt", header=T, https://urldefense.proofpoint.com/v2/url?u=http-3A__as.is&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=O_BK_nzMn3aoOgKF9Hm5nHcz7EiQitAVpsLNjBlwgzs&e=  = T)
row.names(RNAseq) <- RNAseq[,1]
RNAseq <- RNAseq[,-c(1)]
RNAseq <- data.matrix(RNAseq)

#quantile normalization by Dave Tang, source: https://urldefense.proofpoint.com/v2/url?u=https-3A__davetang.org_muse_2014_07_07_quantile-2Dnormalisation-2Din-2Dr_&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=0kYakmXsG_91I29y3gdYAXO6xpzH9gStvGXALbEZzrY&e= 
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#Log 2 transform, this is common for RNA seq analysis:
RNA_log2T <- (log2(RNAseq+1))

#performing quantile normalization, this is not common for RNA seq analysis, but may help for our analysis:
RNA_normal <- quantile_normalisation(RNA_log2T)

#plotting distributions to show the log2 transform, and normalized distributions:

par(mfrow=c(2,1))
boxplot(RNA_log2T, main="After log2 Tranformation")
boxplot(RNA_normal, main= "After log2 Transform & Quantile Normalization")


#separating cases vs control
GE.control <- RNA_normal[,c(seq(1,16))]
https://urldefense.proofpoint.com/v2/url?u=http-3A__GE.case&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=SQ5bRIcImHZJoMdCdfONRZmgYNna3l9WjvJ7OMUkrIQ&e=  <- RNA_normal[,c(seq(17,31))]

#seperating sub-phenotypes (white matter vs grey matter samples)
GE.control.WM <- GE.control[,c(2,4,6,8,10,11,12,13,14,15,16)]
https://urldefense.proofpoint.com/v2/url?u=http-3A__GE.control.GM&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=Vn0OuJBfmWLTeyxriwQPzWZAdj3t1xKoJ0Aiu7IctPo&e=  <- GE.control[,-c(2,4,6,8,10,11,12,13,14,15,16)]

GE.case.WM <- GE.case[,c(2,4,6,8,10,11,12,13,14,15)]
https://urldefense.proofpoint.com/v2/url?u=http-3A__GE.case.GM&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=MFG5E0226a2B6fH3pd8zB-L-pVczFmYdMqYTcCK1q1U&e=  <- GE.case[,-c(2,4,6,8,10,11,12,13,14,15)]


#Only getting the MS-PRGs from Bayesian framework (Andi's work)
MS_PRGs <- read.table("MS_Risk_GeneList_methy.txt", https://urldefense.proofpoint.com/v2/url?u=http-3A__as.is&d=DwIGAg&c=bKRySV-ouEg_AT-w2QWsTdd9X__KYh9Eq2fdmQDVZgw&r=0_Ey8ObnWR6RQ9TmFjrMS3UcaG_BibJ1VjW_J_E1t10&m=HvT-jxUKm5XwhenTSUiJK6EXJFG6QxqRyB1KepUeEek&s=O_BK_nzMn3aoOgKF9Hm5nHcz7EiQitAVpsLNjBlwgzs&e=  = T, header = T)
MS_PRGs <- MS_PRGs$x

ExMatrix <- cbind(GE.case.WM,GE.control.WM)
idx <- which(rownames(ExMatrix) %in% MS_PRGs)
ExMatrix <- ExMatrix[idx,]
sampleLabel <- colnames(ExMatrix)
labels <- c(rep("MS",10), rep("control",11))
names(sampleLabel) <- labels
sampleLabel <- factor(sampleLabel)

#### Running Cogena
devtools::load_all("./")
# KEGG Pathway gene set
annoGMT <- "c2.cp.kegg.v7.01.symbols.gmt.xz"
# GO biological process gene set
# annoGMT <- "c5.bp.v7.0.symbols.gmt.xz"
annofile <- system.file("extdata", annoGMT, package="cogena")

# the number of clusters. It can be a vector.
# nClust <- 2:20
nClust <- 10
# Making factor "Psoriasis" behind factor "ct" means Psoriasis Vs Control
# is up-regualted
sampleLabel <- factor(sampleLabel, levels=c("ct", "Psoriasis"))

# the number of cores.
# ncore <- 8
ncore <- 2

# the clustering methods
# clMethods <- c("hierarchical","kmeans","diana","fanny","som","model",
# "sota","pam","clara","agnes") # All the methods can be used together.
clMethods <- c("hierarchical","pam")


# the distance metric
metric <- "correlation"

# the agglomeration method used for hierarchical clustering
# (hierarchical and agnes)
method <- "complete"

#### Pathway Analysis

# Co-expression Analysis
genecl_result <- coExp(ExMatrix, nClust=nClust, clMethods=clMethods, 
                       metric=metric, method=method, ncore=ncore)

# Enrichment (Pathway) analysis for the co-expressed genes
clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
