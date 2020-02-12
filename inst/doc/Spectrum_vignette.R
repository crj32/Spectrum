## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test1 <- Spectrum(blobs,showpca=TRUE,fontsize=8,dotsize=2)

## ------------------------------------------------------------------------
names(test1)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
RNAseq <- brain[[1]]
test2 <- Spectrum(RNAseq,fontsize=8,dotsize=2)
pca(test2$similarity_matrix,labels=test2$assignments,axistextsize=8,
           legendtextsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
RNAseq <- brain[[1]]
test3 <- Spectrum(RNAseq,showres=FALSE,runrange=TRUE,krangemax=10)

## ------------------------------------------------------------------------
head(test3[[2]]$assignments)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test4 <- Spectrum(brain,fontsize=8,dotsize=2)
kernel_pca(test4$similarity_matrix,labels=test4$assignments,
           axistextsize=8,legendtextsize=8,dotsize=1.5)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
brain1 <- brain[[1]]
brain2 <- brain[[2]]
brain3 <- brain[[3]]
brain1 <- brain1[,-5:-10]
brain_m <- list(brain1,brain2,brain3)
test4 <- Spectrum(brain_m,missing=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test5 <- Spectrum(circles,showpca=TRUE,method=2,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test6 <- Spectrum(spirals,showpca=TRUE,method=2,tunekernel=TRUE,fontsize=8,dotsize=2)

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
test7 <- Spectrum(blobs,FASP=TRUE,FASPk=300,fontsize=8,dotsize=2)

## ------------------------------------------------------------------------
names(test7)

## ------------------------------------------------------------------------
head(test7[[1]])

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
s <- sigma_finder(blobs)
s1 <- ng_kernel(blobs,sigma=s)
e1 <- estimate_k(s1,showplots=FALSE)
r <- cluster_similarity(s1,k=8,clusteralg='GMM')

## ----fig.width=4.5,fig.height=3------------------------------------------
library(Spectrum)
s1 <- CNN_kernel(blobs)
s2 <- CNN_kernel(blobs)
klist <- list(s1,s2)
x <- integrate_similarity_matrices(klist)
e1 <- estimate_k(x,showplots=FALSE)
r <- cluster_similarity(x,k=8,clusteralg='GMM')

## ------------------------------------------------------------------------
## 1. run my clustering algorithm yielding assignments in vector, e.g. 1,2,2,1,2,2...
## 2. reorder data according to assignments
ind <- sort(as.vector(test2$assignments),index.return=TRUE)
datax <- RNAseq[,ind$ix] ## order the original data 
#annonx <- meta[ind$ix,] ## order the meta data
#annonx$cluster <- ind$x ## add the cluster to the meta data
## 3. do heatmap 
# insert your favourite heatmap function

