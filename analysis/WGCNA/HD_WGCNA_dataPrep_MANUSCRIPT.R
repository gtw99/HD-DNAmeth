### load libraries

library(dplyr)
library(stringr)
library(corrplot)
library(cowplot)
library(rgl)
library(WGCNA,options(rlib_downstream_check = FALSE))
library(parallel)

### load normalised betas 
load("Normalised_betas.rdat") 
### load pheno file
pheno <-read.csv("pheno.csv", header = T, stringsAsFactors = F, row.names = 1)

unique(colnames(betas) == pheno$Basename) # sanity check pheno and beta matrix match
identical((pheno$Basename),colnames(betas))

###remove EPIC sex probes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

epicManifest<-read.csv("MethylationEPIC_v-1-0_B4.csv", skip = 7)
rownames(epicManifest)<-epicManifest[,1]
sex_probes <- rownames(epicManifest[which(epicManifest$CHR %in% c("X","Y")),])
betas <- betas[-which(rownames(betas) %in% sex_probes),]


######### Regress out covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Regress out covariates of so the same as PAPER EWAS - 
## STR and EC = Sex, Age, Cell, Plate, Cambridge, London, Manchester, Oxford
## CER = Sex, Age, Plate, Cambridge, London, Manchester, Oxford


# Function to extract residuals

resid <- function(x,Sex, Age, Cell, Plate, Cambridge, London, Manchester, Oxford){
  fit <-try(lm(as.numeric(x) ~ Sex + Age + Cell + Plate + Cambridge + London + Manchester + Oxford), silent=T)
  if(inherits(fit,'try-error')) return(rep(NA,length(sex)))
  return(fit$residuals + fit$coefficients[1])
} 


Sex<-as.factor(pheno$Sex)
Age <- as.numeric(pheno$Age)
Cell <- as.numeric(pheno$prop)
Plate <- as.numeric(pheno$Plate)
Cambridge <- as.factor(pheno$Institute_Cambridge)
London <- as.factor(pheno$Institute_London)
Manchester <- as.factor(pheno$Institute_Manchester)
Oxford <- as.factor(pheno$Institute_Oxford)


cl<-makeCluster(16)
normMat <- t(parApply(cl,betas,1,resid, Sex, Age, Cell, Plate, Cambridge, London, Manchester, Oxford))
stopCluster(cl)
colnames(normMat) <- colnames(betas)
rownames(normMat) <- rownames(betas) 
norm <- normMat



# Filter to most variable residualised probes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# perform individaully on each brain region


varTest <- function(x){
  mad <- mad(x)   ### median absolute deviation
  med <- median(x)
  return(c(mad,med))
}

cl<- makeCluster(16)
resVar <-t(parApply(cl, norm, 1, varTest))
stopCluster(cl)
median(resVar[,1])

med <- rownames(norm[which(resVar[,1] > median(resVar[,1])),])

### IDs of all the variable probes from each region
allProbes <- unique(c(med,med_region,med_region))

### removed non nariable probes (probes with variance > the region with lowest median)

t(norm[allProbes,])->dat


#### Detect and remove outliers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

options(stringsAsFactors = FALSE)

### Check for extreme outliers 
### Cluster samples on Euclidean distance
sampleTree = hclust(dist(dat), method = "average");

### Clustering dendrogram:

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
abline(h = cutHeight_val, col = "red") # Can set arbitrary cutoff threshold to observe which samples may be outliers 

### generate principal components from variable betas 
pc <- prcomp(dat)

### visual inspection of outliers using PCs
par(mfrow=c(2,2),mar = c(4,4,4,4))
plot(pc$x[,1],pc$x[,2],col = as.factor(rownames(dat) %in% rownames(dat)),pch = 4,xlab = "PC1", ylab = "PC2")
legend("bottomleft", c("Exclude", "Keep"), col = c("black", "red"), pch = 1)
plot(pc$x[,2],pc$x[,3],col = as.factor(rownames(dat) %in% rownames(dat)),pch = 4,xlab = "PC2", ylab = "PC3")
plot(pc$x[,1],pc$x[,4],col = as.factor(rownames(dat) %in% rownames(dat)),pch = 4,xlab = "PC1", ylab = "PC4")
plot(pc$x[,3],pc$x[,1],col = as.factor(rownames(dat) %in% rownames(dat)),pch = 4,xlab = "PC3", ylab = "PC1")


### function to identify samples to keep
cutNkeep <- function(sampleTree, cutnum, region, data){
  clust <- cutreeStatic(sampleTree, cutHeight = cutnum, minSize = 10)
  message(paste(table(clust)))
  keepSamples = (clust==1)
  keepFrame <- data[keepSamples,]
  message(paste(region,"has dimensions", dim(keepFrame)))
  assign(paste("dat",region,sep = ""), keepFrame, envir = parent.frame())
  assign(paste("keepSamples",region,sep = ""), keepSamples, envir = parent.frame() )
}

cutNkeep(sampleTree,cutHeight_val,dat)

### keep only non-outlying samples 
clust = cutreeStatic(sampleTreeStri, cutHeight = cutHeight_val, minSize = 10)
table(clust)
keepSamples = (clust==1)
dat = dat[keepSamples, ]
dim(dat)

