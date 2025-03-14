### load libraries 

library(plyr) 
library(dplyr) 
library(stringr)
library(corrplot)
library(cowplot)
library(rgl) 
library(WGCNA)  
library(parallel)
library(lm.beta)   
library(utils) 
library(tibble)
library(clusterProfiler)
library(KEGGREST)
library(missMethyl) 
library(rrvgo)



load("blockwiseMods.Rdata")
### load dat file containing only the most variable probes and with outlying samples removed
load(file="dat.RData")
pheno<-read.csv("pheno.csv", header = T, stringsAsFactors = F, row.names = 1)



moduleColors = labels2colors(bwnet.$colors)
table(moduleColors)
dg <- which(moduleColors == "col") #### take module of interest if using whole module

### if examining hub probes run these line instead:
# MM_col <- read.csv("MM_col_file.csv", stringsAsFactors = FALSE, row.names = 1)
# hub_probes <- subset(MM_col abs(moduleMembership) > 0.8 & geneP < 0.05)
# dg <- row.names(hub_probes)

dg_names<- names(bwnet$colors[dg])
all_names <- names(bwnet$colors)


### run KEGG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gstKEGG <- gometh(sig.cpg=dg_names, all.cpg=all_names, collection="KEGG", array.type = "EPIC",
                  plot.bias=TRUE)
topRcolKEGG <- topGSA(gstKEGG, n=1000)

### plot KEGG

goSet <- topRcolKEGG 
goSet$Description  <- factor(goSet$Description , levels = goSet$Description )

ggplot(goSet[1:2,], aes(x = as.factor(Description), y = -log10(P.DE),size = DE/N))+
  geom_point(fill = "lavenderblush3",colour = "black", shape = 21)+
  coord_flip()+
   ylab("-log10(P)")+
   xlab("KEGG Term")+
labs(size="Proportion")+
theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.y = element_text(size = 12)) +
  theme_bw()


### run GO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gst <- gometh(sig.cpg=dg_names, all.cpg=all_names, collection="GO", array.type = "EPIC",
              plot.bias=TRUE)
topRcol <- topGSA(gst, n=1000)

### plot GO

goSet <- topRcol
goSet$TERM <- factor(goSet$TERM, levels = goSet$TERM)

ggplot(goSet[1:10,], aes(x = as.factor(TERM), y = -log10(P.DE),size = DE/N))+
  geom_point(fill = "lavenderblush3",colour = "black", shape = 21)+
  coord_flip()+
  ylab("-log10(P)")+
   xlab("Go Term")+
labs(size="Proportion")+
theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.y = element_text(size = 12)) +
  theme_bw()



#### collapse GO terms with rrvgo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gst$ID <- rownames(gst) 

scores <- setNames(-log10(gst$P.DE), gst$ID)
simMatrixBP <- calculateSimMatrix(rownames(gst),
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
reducedTermsBP <- reduceSimMatrix(simMatrixBP,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
simMatrix <- calculateSimMatrix(rownames(gst),
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")
reducedTermsCC <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
simMatrix <- calculateSimMatrix(rownames(gst),
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")
reducedTermsMF <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")

reducedTerms <- rbind(reducedTermsCC,reducedTermsMF,reducedTermsBP)

gst$parentTerm <- gst$TERM
gst[rownames(reducedTerms),"parentTerm"] <- reducedTerms$parentTerm
gst$Sort <- as.numeric(!gst$TERM == gst$parentTerm)
gst <- gst[order(gst$parentTerm,gst$Sort),]
gst$plotIndex <- 1:nrow(gst)

gst <- gst[which(gst$TERM%in% unique(gst$parentTerm)),]
topRcol <- topGSA(gst, n=1000)

### plot reducded GO

goSet <- topRcol
goSet$TERM <- factor(goSet$TERM, levels = goSet$TERM)

ggplot(goSet[1:10,], aes(x = as.factor(TERM), y = -log10(P.DE),size = DE/N))+
  geom_point(fill = "lavenderblush3",colour = "black", shape = 21)+
  coord_flip()+
  ylab("-log10(P)")+
   xlab("Go Term")+
labs(size="Proportion")+
theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.y = element_text(size = 12)) +
  theme_bw()
