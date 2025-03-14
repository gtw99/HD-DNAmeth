library(WGCNA)
library("tidyverse")
library(plyr) 
library(dplyr) 
library(stringr)
library(corrplot)
library(cowplot)  
library(parallel)
library(lm.beta)   
library(utils) 
library(tibble)
library(reshape2)
library(gridExtra)
library("ggplotify")

### load dat file containing only the most variable probes and with outlying samples removed
load(file="dat.RData")
### load pheno file for trait association
pheno <-read.csv("pheno.csv", header = T, stringsAsFactors = F, row.names = 1)


### Identify the soft threshold power ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

enableWGCNAThreads(8)
powers = 1:20
### Call the network topology analysis function
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, networkType = "unsigned")

### plot results

## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"),ylim=c(0.2,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.8, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
abline(h = 0.8, col = "orange")
##Mean connectivity as
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.8, col = "red")


### Generate the modules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

softPower = value; 

enableWGCNAThreads(32)
bwnet = blockwiseModules(dat, maxBlockSize = 10000, power = softPower, TOMType = "unsigned", deepSplit = 0,
 minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3)
save(bwnet, file = "blockwiseMods.Rdata")


### Trait association ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### match pheno to samples in dat (i.e. remove outliers from pheno) and reduce to traits to test
pheno <- pheno[match(rownames(dat), pheno_Striatum$Basename),c("Sample_ID", "Basename", "Phenotype", "Sex", "Age","prop","Institute_London", "Institute_Manchester", "Institute_Cambridge", "Institute_Oxford")]

pheno<- pheno%>% 
       rename("NeuN" = "prop",
       "LNDBB" = "Institute_London",
       "MBB" = "Institute_Manchester",
       "CBB" = "Institute_Cambridge",
       "OBB" = "Institute_Oxford")


####### Function to perform module trait associations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# tVec: Vector of labels for trait assocation, must be in form binary or quantitative
# tFile: Data frame with tVec labels as colnames (pheno file)
# mFile: module file containing module eigengenes
# rDat: "Raw" base gene expression / methylation matrix used to construct module

traitModCor <- function(tVec,tFile,mFile,rDat){
  message(paste("You've started the function, at least"))
  # Define your variables
  nGenes = ncol(rDat)
  nSamples= nrow(rDat)
  moduleColors = labels2colors(mFile$colors)
  MEs0= moduleEigengenes(rDat, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  message(paste("You've gotten to the pre-check"))
  # 
  if(!identical(rownames(rDat),rownames(MEs))){
    
    stop("Raw data matrix and module eigengenes are discordant")
    
  } 
  
  if(length(unique(tVec %in% colnames(tFile))) > 1 | unique(tVec %in% colnames(tFile)) == FALSE){
    
    stop("Label vector and trait file are discordant")
    
  }
  
  message(paste("You've gotten to the results file creation"))
  
  moduleTraitCor <- matrix(data = NA, ncol = length(tVec), nrow = ncol(MEs))
  colnames(moduleTraitCor)<- tVec
  rownames(moduleTraitCor) <- colnames(MEs)
  
  moduleTraitPvalue <- matrix(data = NA, ncol = length(tVec), nrow = ncol(MEs))
  colnames(moduleTraitPvalue)<-  tVec
  rownames(moduleTraitPvalue) <- colnames(MEs)

  message(paste("You've gotten to the correlation loop"))
  
  for(t in tVec){
    if(length(unique(tFile[,t])) > 2){ # if Binary run spearman cor
      for(m in colnames(MEs)){
        
        try(res <-cor.test(MEs[,m] , as.numeric(tFile[,t]),method = "pearson"), silent = TRUE)
        if(class(res) != "try-error")
          moduleTraitCor[m,t]<-as.numeric(res$estimate)
        moduleTraitPvalue[m,t]<-res$p.value
      }
      
    } else { # do the same as above but if quantitative then run different pearson cor
      for(m in colnames(MEs)){
        try(res <-cor.test(MEs[,m] , as.numeric(tFile[,t]),method = "spearman"), silent = TRUE)
        if(class(res) != "try-error")
          moduleTraitCor[m,t]<-as.numeric(res$estimate)
        moduleTraitPvalue[m,t]<-res$p.value
      } 
    }
    
    message(paste("Testing",t))
  }
  assign("moduleTraitCor", moduleTraitCor, envir = parent.frame() )
  assign("moduleTraitPvalue", moduleTraitPvalue, envir = parent.frame() )
}


names(bwnet)
moduleColors = labels2colors(bwnet$colors)

labs <- colnames(pheno[c(3:10)])  


traitModCor(tVec = labs, tFile = pheno, mFile = bwnet, rDat = dat) 


### remove confounded modules (those associated with confounding variables)
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
rem <- which(moduleTraitPvalue$Sex < 0.05 | moduleTraitPvalue$Age < 0.05 | moduleTraitPvalue$NeuN < 0.05| 
  moduleTraitPvalue$LNDBB < 0.05 | moduleTraitPvalue$MBB < 0.05 | moduleTraitPvalue$CBB < 0.05 | moduleTraitPvalue$OBB < 0.05 |  rownames(moduleTraitPvalue) == "MEgrey")
#Remove confounded mocules
moduleTraitPvalue <- moduleTraitPvalue[-rem,]
moduleTraitCor<- moduleTraitCor[-rem,]


moduleTraitPvalue <- as.matrix(moduleTraitPvalue)

### create Heatmap to display correlations and their p-values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3,6.5,3,1));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = labs,
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))          


### t.test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### get the modules that are not assocaited with coVars
modKeep <- rownames(moduleTraitPvalue)

MEs0= moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
dim(MEs)

t(MEs)->MET
moduleCounts <- as.data.frame(table(moduleColors))


###  t test function  

  ttest <- function( row, Diag){     
    
    ttest.result <- t.test( row ~ factor(Diag) )
    
    
    
    return(ttest.result$p.value)
  }

  
Diag     <- factor(pheno$Phenotype) 
 
diagtest <-  t(apply( MET, 1,ttest, Diag))  ### apply t-test function
diagtest <- t(diagtest)


###filter out modules associated with coVars
diagtest <- as.data.frame(diagtest[row.names(diagtest) %in% modKeep, ], row.names = modKeep,)
diagtest <- as.matrix(diagtest)
pOrder <- order(diagtest[,1])
diagtest2 <- as.data.frame(diagtest[pOrder,])
head(diagtest2)

colnames(diagtest2)[colnames(diagtest2) == 'diagtest[pOrder, ]'] <- 'pValues'
as.data.frame(MET)->MET
as.matrix(MET['ME_selected_moule_colour',])-> ME_module
t(ME_module)->ME_module

### box plot of module eigengenes
ggplot(pheno, aes(x = as.factor(Phenotype), y = ME_module, colour = as.factor(Phenotype),fill = as.factor(Phenotype)))+
geom_boxplot(color = "black",lwd = 0.8, outlier.shape = NA)+
geom_jitter(width = 0.3, shape = 21, colour = "black", fill = "black")+
theme_bw()+
ylab("Module eigengene")+
xlab(NULL)+
ggtitle("red")+ 
theme(legend.position = "none")+
theme(plot.title = element_text(hjust = 0.5, size = 18),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14))+
scale_x_discrete(labels=c("1" = "Control",
                              "2" = "HD"))

### module memebership  analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### change bwnet.XXX as necessary
##### change 'dat' as necesssary 
##### change pheno as necessary 
##### change module of interest as necessary

moduleColors = labels2colors(bwnet$colors)
table(moduleColors)

MEs0= moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

rownames(MEs) == rownames(dat) 

weight = as.data.frame(as.numeric(pheno$Phenotype));
names(weight) = "Phenotype"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(pheno)));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(dat, weight, use = "p",method = "spearman"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(pheno)));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
head(MMPvalue)

dg <- which(moduleColors == "MEcol") #### take module of interest 
dg_names<- names(bwnet$colors[dg])
all_names <- names(bwnet$colors)

epicManifest<-read.csv("MethylationEPIC_v-1-0_B4.csv", skip = 7)
rownames(epicManifest)<-epicManifest[,1]

plotMM <- data.frame(row.names = rownames(GSPvalue),moduleMembership = geneModuleMembership$MMcol,MMPvalue = MMPvalue$p.MMcol,
                     geneCor = geneTraitSignificance$GS.Phenotype,geneP = GSPvalue$p.GS.Phenotype,epicManifest[colnames(dat[,all_names]),])


### make summary with rank for indivdual modul then add with EPIC annotation

plotMM <- plotMM[dg_names,] ### cpg names in the module dg names is first generated above for the trait associatons 

#### add in MMrank and GSrank before EPIC annotation
plotMM <- plotMM %>%
  add_column(MMrank = rank(abs(plotMM$moduleMembership)),
              GSrank = rank(-log10(plotMM$geneP)),
             .before = "IlmnID")
             
### summative rank     
sumrank = plotMM$MMrank + plotMM$GSrank
plotMM <- plotMM %>%
  add_column(sumrank = sumrank,
             .before = "IlmnID")

dim(plotMM)
### run correlation test
est <- cor.test(abs(plotMM$moduleMembership),-log10(plotMM$geneP))$estimate ### negative log ten
pval <- cor.test(abs(plotMM$moduleMembership),-log10(plotMM$geneP))$p.value
MMPScorrelation<-data.frame(est,pval)

### plot ### 
ggplot(plotMM, aes(x = abs(moduleMembership), y = -log10(geneP), col = sumrank))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = 0.8)+
  geom_smooth(method = "lm")+
  theme_bw()+
  ylab("Probe Significance for HD (-log10(p))")+
  xlab("Module Membership")+
  ggtitle("lavenderblush3", subtitle = "corr = 0.458, p =1.91e-12")+ 
  theme(plot.title = element_text( size = 22),
  plot.subtitle = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16))+
  theme(legend.position = c(0.1, 0.85),
  legend.text = element_text(size = 14),
   legend.title = element_text(size = 14))+
  scale_color_gradient(low="palevioletred2", high="palevioletred4")

#### Hub genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

top_MM <- plotMM[order(plotMM[, 7],decreasing = TRUE), ]   ##### orders membership summary file by MM

### select hub probes 
hub_probes <- subset(top_MM, abs(moduleMembership) > 0.8 & geneP < 0.05) 

testDF <- hub_probes
testDF$UCSC_RefGene_Name <- as.character(testDF$UCSC_RefGene_Name)

#### seperate out UCSC multiple gene names into seprate columns

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(rownames(testDF)), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")

for(i in 1:length(rownames(testDF))){
  cpgIndex <- rownames(testDF)[i]
  test1 <- as.character(testDF[cpgIndex,"UCSC_RefGene_Name"])
 
if(str_detect(test1,";") == TRUE){
  strIndex <- as.data.frame(str_locate_all(test1,";"))$end
  length(strIndex)
  strIndex <- c(0,strIndex,str_length(test1))
  storage <- c()
  for(y in strIndex){
    if(which(strIndex == y) == 2){
      storage <- append(storage,str_sub(test1,
                                 start = strIndex[which(strIndex == y) -1],
                                 end = y - 1))
    }else if(y >0 & y != max(strIndex)){
      storage <-append(storage,str_sub(test1,
                                 start = strIndex[which(strIndex == y) -1]+ 1,
                                 end = y - 1))
    }else if(y == max(strIndex)){
      storage <-append(storage,str_sub(test1,
                                       start = strIndex[which(strIndex == y) -1]+ 1))
    } 
  }
    sumGene <- as.character(unique(storage))
    
    if(length(sumGene) == 1){
      testDF[cpgIndex,"UCSC1"] <- sumGene
    }else if(length(sumGene) == 2){
      testDF[cpgIndex,c("UCSC1","UCSC2")] <- sumGene[1:2]
    }else if(length(sumGene) > 2){
      testDF[cpgIndex,c("UCSC1","UCSC2","UCSC3")] <- sumGene[1:3]
    }
  } else if(str_detect(test1,";") == FALSE){
  testDF[cpgIndex,"UCSC1"] <- as.character(testDF[cpgIndex,"UCSC_RefGene_Name"])
  } 
setTxtProgressBar(pb, i)
}
                               
write.table(testDF, file = "hub_probes.txt")

probes_per_gene <-table(testDF$UCSC1)
probes_per_gene_sort1 <- probes_per_gene[order(probes_per_gene)] 
write.table(probes_per_gene_sort1, file = "hub_probes_per_gene.txt")