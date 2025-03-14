### load libraries ~~~~~~~~~~~~~~~~~~~~

library(corrplot)
library(Hmisc)
library(wateRmelon)
library(methylumi)
library(plyr)#for relvelling
library(parallel) #for cluster
library(Haplin) #for pQQ
library(qqman) #for manhattan
library(dplyr)
library(bacon)
library(ggrepel)
library(stringr)
library(tibble)


### load methylation betas and pheno file 
load("betas.rdat") 
pheno <-read.csv("pheno_file.csv", header = T, stringsAsFactors = F, row.names = 1)


### Run the linear model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### LM function
### For STR and EC the cell proportion is included as a covariate. For the CER it is removed.

lm_fun <- function(betas,Phenotype, Sex, Age, Cell, Plate, Cambridge, London, Manchester, Oxford){
  res<-lm(betas ~ Phenotype + Sex + Age + Cell + Plate + Cambridge + London + Manchester + Oxford)
  return(cbind(coef(summary(res))[2,]))
}

cl<- makeCluster(16)

Phenotype<-pheno$Phenotype
Sex<-pheno$Sex
Age <- pheno$Age
Cell <- pheno$prop
Plate <- pheno$Plate
Cambridge <- pheno$Institute_Cambridge
London <- pheno$Institute_London
Manchester <- pheno$Institute_Manchester
Oxford <- pheno$Institute_Oxford

model<-t(parApply(cl,betas,1,lm_fun, Phenotype, Sex, Age, Cell, Plate, Cambridge, London, Manchester, Oxford))

stopCluster(cl)

colnames(model)<-c("Est","SE","Tval","Pval") #rename blank columns to summary statistics
head(model)
model <- as.data.frame(model)


### Calculating lambda ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

lambda1 <- LambdaInf(na.omit(model$Pval))

### bacon correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### for brain regions with lambda > 1.2

#### need effect size and standard errors for bacon
es<-modelSTR[,"Est"]
se<-modelSTR[,"SE"]
identical(rownames(es), rownames(se))
library(BiocParallel)
register(MulticoreParam(32, log=TRUE))

bc <- bacon(NULL, es, se)
es.b<-es(bc)
se.b<-se(bc)


estimates(bc)
write.table(estimates(bc), file = "model_bacon2_estimates", row.names = F)

plot(bc, type="hist")
lambdaBacon <- LambdaInf(na.omit(pval(bc)))
pval.bc <- pval(bc)
plot(bc, type="qq")


### annotatation and plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

epicManifest<-read.csv("MethylationEPIC_v-1-0_B4.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(model), epicManifest$Name),]

### bind bacon results and desired EPIC manifest metrics
bac_anno<-cbind(pval.bc,es.b,se.b, epicManifest_STR[,c("EPIC annotation columns", "...", "...")])

model <- cbind(model, bac_anno)

#Recode x/y
model$CHR<-as.character(model$CHR)
model$CHR[which(model$CHR == "X")]<-23
model$CHR[which(model$CHR == "Y")]<-24
model$CHR<-as.numeric(model$CHR)
model<-model[which(model$CHR != ""),]### to allow manhatten plot

model  <- model [order(model$pval.bc), ]   ### reorder for USCS funciton below
head(model)

#### seperates out UCSC columns currently just the top 50 to allow annotation on plots
testDF <- model

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(rownames(testDF)), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")

for(i in 1:50){
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

#### replaces gaps and NAs with cpgnames
testDF[which(testDF$UCSC1 == ""),"UCSC1"] <- rownames(testDF[which(testDF$UCSC1 == ""),])
testDF[which(is.na(testDF$UCSC1)),"UCSC1"] <- rownames(testDF[which(is.na(testDF$UCSC1)),])


###need dataframe in this format for ggplot
RESULTS <- testDF
colnames(RESULTS)[c(5,10)] <- c("P","BP")  ### need to ensure right columns selected for P and BP -
RESULTS <- RESULTS[,c("CHR","BP","P","UCSC1")] 
RESULTS$CHR <- as.numeric(RESULTS$CHR)
RESULTS$BP <- as.numeric(RESULTS$BP)
RESULTS <- na.omit(RESULTS)

don <- RESULTS %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(RESULTS, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>% 
  mutate( is_annotate=ifelse(P < 6.271511e-08, "yes", "no"))
  
#rownames(don) <- don$SNP
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

### Manhatten plot 
ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("powderblue", "steelblue3"), 22 )) +
    geom_hline(yintercept = -log10(6.271511e-08), linetype = 'dotted', col = 'red')+
    geom_hline(yintercept = -log10(1.0e-05), linetype = 'dotted', col = 'black')+
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,11) ) +     # remove space between plot area and x axis
    # Custom the theme:
    theme_bw() +
  
    # Add label using ggrepel to avoid overlapping, remove repel as necessary - label adds box, text just adds text
    geom_text_repel( data=subset(don, is_annotate=="yes"), aes(label=UCSC1), size=2)+
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y=element_text(size=rel(0.5)),
      axis.text.x=element_text(size=rel(0.5), angle = 45)
    )+
	xlab("Chromosome")

### Effect size comparisons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Brain regions

load('EC_model.RData')
load('Cer_model.RData')


modelSTR_100 <- model[1:100,]

# # # EC # # #

### pull out the 100 most significant STR probes from the EC data
STR100_inEC <-modelEC[match(rownames(modelSTR_100), rownames(modelEC)),]

### run correltaion
EC_ES_corr<-rcorr(modelSTR_100$es.b, STR100_inEC$es.b, type="pearson")
corr_str_EC <-table(EC_ES_corr)

top_stri_EC<-cbind(modelSTR_100$es.b, STR100_inEC$es.b)
colnames(top_stri_EC)[c(1,2)] <- c("ES of STR CpGs","ES of the STR CpGs in the EC")
top_stri_EC <- as.data.frame(top_stri_EC)
plot(top_stri_EC, xlab = "ES of STR CpGs", ylab = "ES of the STR CpGs in the EC")
abline(v = 0, col="red", lty=5)
abline(h = 0, col="red", lty=5)

### run binominal test
bires_str_EC <- binom.test(nrow(modelSTR_100[modelSTR_100$es.b > 0 & STR100_inEC$es.b > 0,]) + 
nrow(modelSTR_100[modelSTR_100$es.b < 0 & STR100_inEC$es.b< 0,]), 100)


# # # CER # # # 
### pull out the 100 most significant STR probes from the CER data
STR100_inCer <-modelCer[match(rownames(modelSTR_100), rownames(modelCer)),]

### run correlation
Cer_ES_corr<-rcorr(modelSTR_100$es.b, STR100_inCer$Est, type="pearson")
corr_str_Cer <-table(Cer_ES_corr)

top_stri_Cer<-cbind(modelSTR_100$es.b, STR100_inCer$Est)
colnames(top_stri_Cer)[c(1,2)] <- c("Striatum top 100 ES","Striatum top 100 probes in the Cer ES")
top_stri_Cer <- as.data.frame(top_stri_Cer)
plot(top_stri_Cer, xlab = "ES of STR CpGs", ylab = "ES of the STR CpGs in the CER")
abline(v = 0, col="red", lty=5)
abline(h = 0, col="red", lty=5)

### run binominal test
bires_str_Cer <- binom.test(nrow(modelSTR_100[modelSTR_100$es.b > 0 & STR100_inCer$Est > 0,]) + 
nrow(modelSTR_100[modelSTR_100$es.b < 0 & STR100_inCer$Est< 0,]), 100)



### Comparison with horvath et al 2016

### load horvath data 
horvath <- read.csv("Horvath_20165_HD_EWAS_summary_stats.csv", header = T, stringsAsFactors = F, row.names = 1)

horvath_order <- horvath[order(horvath$p.metaAnalysis), ]

#subset common probes 
hor_com_100  <- subset(horvath, rownames(horvath) %in% rownames(modelSTR_100)) 
STR_com_100  <- subset(modelSTR_100, rownames(modelSTR_100) %in% rownames(horvath))
hor_com_100 <-hor_com_100 [match(rownames(STR_com_100), rownames(hor_com_100)),] ## reorders the horvath dataset 
identical(rownames(STR_com_100), (rownames(hor_com_100)))

STR_com_100$es.b <- 100*(STR_com_100$es.b) ###  ES need to be as percentage to match Horvath

### run correlation 
ES_corr <-rcorr(STR_com_100$es.b, hor_com_100$Z.metaAnalysis, type="pearson")
corr_ES <-table(ES_corr)

plot(STR_com_100$es.b, hor_com_100$Z.metaAnalysis, xlab = "STR Probe ES" , ylab = "Horvath et al. Z-scores")
abline(v = 0, col="red", lty=5)
abline(h = 0, col="red", lty=5)


### run binominal test
bi_STR_hor <- binom.test(nrow(STR_com_100[STR_com_100$es.b > 0 & hor_com_100$Z.metaAnalysis > 0,]) + 
nrow(STR_com_100[STR_com_100$es.b < 0 & hor_com_100$Z.metaAnalysis< 0,]), 48)


