###Load packages 
library(methylumi)
library(wateRmelon)
require(gdata)
library(minfi)
library(ggplot2)
library(gdata)
require(gridExtra)
require(IlluminaHumanMethylationEPICmanifest) 
###EPIC manifest which includes info on probe location, probe type etc. 
library(tidyr)
library(dplyr)
library(cets)


### Read in pheno file and assign
pheno<-read.csv("pheno_file")

### sort sample names
pheno$Basename <- as.character(pheno$Basename) ### ensures chip name is character vector and assigns
pheno$Basename2<-pheno$Basename
pheno<-separate(data = pheno, col = Basename2, into = c("SentrixID", "Position"), sep="_")

###assign  idats to a variable              
idatPath<-c("IDAT_folder_path")

### Create methylumiSet and RGset objects
msetEPIC <- readEPIC(idatPath=idatPath, barcodes=pheno$Basename, parallel = FALSE, force=T)
RGset <- read.metharray.exp(base = idatPath, targets = pheno, force = TRUE)

### checks the basename ID matches mset and RGset
identical((pheno$Basename),colnames(msetEPIC))
identical((pheno$Basename),colnames(RGset))


### Created copy of pheno file to bind QC results
QCmetrics<-pheno

###Boolean record of which samples have already failed
SamplesFail<-as.logical(rep("FALSE", nrow(pheno)))
###Entries will be changed to TRUE as samples fail
Stepsummary<-as.data.frame(matrix(ncol=0, nrow=2))
rownames(Stepsummary)<-c("Failed This Step", "Total Failed")


### Check Signal Intensities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###The intensity check is the biggest indicator of sample quality. 

###Generating methylated (M) and unmethylated (U) intensities.
m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)

#Calculate the median of said intensities per sample.
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)
 
QCmetrics<-cbind(pheno,M.median, U.median)

### plot intensities 
### colour by desired variable eg bisulfite plate
plotfactor<-factor(pheno$variable, levels=c(unique(pheno$variable), "FullyMethylated"))
plotfactor[pheno$Control]<-"FullyMethylated"
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity", main="Histogram of Median Methylated Intensities", cex.main=0.7)
hist(U.median, xlab = "Median U intensity", main="Histogram of Median Unmethylated Intensities", cex.main=0.7)
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col = rainbow(nlevels(plotfactor))[factor(plotfactor)])
par(xpd=TRUE)
legend("topleft", levels(factor(plotfactor)), col = rainbow(nlevels(plotfactor)), pch = 10, cex=0.5)

#check that the order of the intesities matches the pheno 
as.factor(colnames(m_intensities))->colnames(m_intensities)
m_intensities<-m_intensities[order((pheno$Basename))]
identical(colnames(m_intensities),pheno$Basename)

as.factor(colnames(u_intensities))->colnames(u_intensities)
u_intensities<-u_intensities[order((pheno$Basename))]
identical(colnames(u_intensities),pheno$Basename)

###Remove fully methylated controls and save intensities 

FMpheno<-QCmetrics[pheno$Control,c("Basename", "SentrixID", "Position", "M.median", "U.median")]
### time stamp fullymeth controsl
times<-data.frame(Date_Ran=c("DATE","DATE"), Time_Ran=c("time","time"))
#study info
info<-data.frame(Study=rep("HD_EWAS",nrow(FMpheno)), iDAT_Location=rep(idatPath, nrow(FMpheno)))
FMpheno<-cbind(info, times, FMpheno)
write.csv(FMpheno, "FullyMethylatedControlSamples.csv", row.names=FALSE)

###remove methylated controls from all variables
M.median<-M.median[!pheno$Control]
U.median<-U.median[!pheno$Control]
msetEPIC<-msetEPIC[,!pheno$Control]
RGset<-RGset[,!pheno$Control]
SamplesFail<-SamplesFail[!pheno$Control]
QCmetrics<-QCmetrics[!pheno$Control,]
pheno<-pheno[!pheno$Control,]

###assign samples below a threshold - 1000 for this study
lowintensitysamples<-which(M.median < 1000 | U.median < 1000)
###plot intensity threshold 
Intensity<-rep("OK", nrow(pheno))
Intensity[lowintensitysamples] <-"LowIntensity"
plotfactor<-as.factor(Intensity)
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col=rainbow(2)[factor(plotfactor)])
abline(v = 1000, col = "red")
abline(h = 1000, col = "red")
legend("topleft", levels(factor(plotfactor)), pch = 16, col = rainbow(2))

SamplesFail[which(Intensity=="LowIntensity")]<-TRUE
QCmetrics<-cbind(QCmetrics, Intensity)
Step1<-c(sum(Intensity=="LowIntensity"),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step1)


###Bisulphite conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### A bisulphite conversion statistic for each sample was calculated, and a histogram of the results plotted.
### Samples with a conversion < 80% fail the QC

Bisulphite<-bscon(msetEPIC)
hist(Bisulphite, xlab = "Median % BS conversion", main = "Histogram of Bisulphite Conversion Statistics")
QCmetrics<-cbind(QCmetrics, Bisulphite)
SamplesFail[which(Bisulphite<80)]<-TRUE
Step2<-c(sum(Bisulphite<80, na.rm=T),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step2)


### predict sex of samples using minifi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GRset <- mapToGenome(RGset) 
PredictedSex_minifi <- getSex(GRset, cutoff = -2) 
Sexplot_data <- as.data.frame(PredictedSex_minifi) 

###Sexplot_data <- Sexplot_data[order(row.names(Sexplot_data)),] 
PredictedSex <-Sexplot_data$predictedSex 
QCmetrics <- cbind(QCmetrics,PredictedSex) 

### add the reported sex data 

Sexplot_data <- cbind(Sexplot_data, QCmetrics$Sex) 
colnames(Sexplot_data)[colnames(Sexplot_data)=="QCmetrics$Sex"] <- "Reported Sex" 

#replace blanks with NA 
Sexplot_data$`Reported Sex`[Sexplot_data$`Reported Sex` == ""] <- "NA" 

ggplot(Sexplot_data, aes(Sexplot_data$xMed,Sexplot_data$yMed,                         
                         colour = Sexplot_data$`Reported Sex`)) +   
  geom_point() + 
  labs(x= "X Chr, median total intensity (log2)", 
       y ="Y Chr, median total intensity (log2)", 
       colour = "Reported Sex") 

ReportedSex <- as.character(QCmetrics$Sex) 
QCmetrics$MismatchSex<-PredictedSex!=ReportedSex 

SamplesFail[which(PredictedSex!=ReportedSex)]<-TRUE 
Step3<-c(length(which(PredictedSex!=ReportedSex)),sum(SamplesFail)) 
Stepsummary<-cbind(Stepsummary,Step3) 


### Genetic correlations id samples from same individual ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### pulls out the unique IDs for each individual and assigns 
brainIDs<-unique(pheno$Individual_ID) # n = 42
###converts Individual IDs to table and counts
singles <- table(as.character(pheno$Individual_ID))
### names individual ids with only 1 count and assigns
singleID <- names(singles[which(singles[] == 1)])
###removes singles 
brainIDs<-setdiff(brainIDs, singleID) # 41
minbrainIDcor<-NULL
numbersamples<-NULL

Individual_IDs<-unique(as.character(pheno$Individual_ID))
Sample_IDs<-unique(pheno$Sample_ID) 
###Use max corr not min as min would cause 3 samples to fail 
###in  a triplicate whereas max ensures two samples with a good correlation pass when the third sample has poor correlations with both

snpCor2 <- snpCor[-which(pheno$Individual_ID==singleID),-which(pheno$Individual_ID==singleID)]
pheno2 <- pheno[-which(pheno$Individual_ID==singleID),]
Individual_IDs<-unique(pheno2$Individual_ID)
Sample_IDs<-unique(pheno2$Sample_ID)

maxrelatedcors <- matrix(data=NA, nrow=1, ncol=2)
colnames(maxrelatedcors) <- c('Basename','Corr')

for(j in 1:length(Individual_IDs)){
  ID <- Individual_IDs[j]
  IDs <- which(QCmetrics$Individual_ID==ID)
  relatedcors <- snpCor[IDs,IDs]
  for(i in 1:length(IDs)){
    maxcor<-max(relatedcors[i,], na.rm = TRUE)
    maxrelatedcors <- rbind(maxrelatedcors, c(rownames(relatedcors)[i], maxcor))
  }
}
maxrelatedcors <- maxrelatedcors[-1,]
maxrelatedcors[,'Corr'] <- as.numeric(maxrelatedcors[,'Corr'])
MaxCor <- maxrelatedcors[match(QCmetrics$Basename,maxrelatedcors[,'Basename']),'Corr']
MaxCor <- (as.numeric(MaxCor))

hist(MaxCor, main="Maximum correlation in samples from the same individuals", xlab="Max Correlation")

print(paste("Samples with low correlation to samples from same individual: ", QCmetrics[which(MaxCor<0.8), 'Sample_ID'], sep=""))
####
MaxCor <- (as.numeric(MaxCor))

QCmetrics<-cbindX(QCmetrics,as.data.frame(MaxCor))
SamplesFail[which(MaxCor< 0.8)]<-TRUE

Step4<-c(sum(QCmetrics$MaxCor < 0.8, na.rm=T), sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step4)


###Hierarchical clustering to confirm relatedness
ind_bn<-substr(rownames(msetEPIC),1,2)
which(ind_bn=="rs")->rsbn
betas(msetEPIC)->beta_COMB
beta_COMB[rsbn,]->rs_bn
dim(rs_bn)
#tellsmefeatures=65andsamples=384.Thisisexpectedas65rsprobes
na.omit(rs_bn)->rs_bn2
cor(rs_bn2)->cornaomit
write.csv(cornaomit,"cornaomit_v2.csv")
#dohistogramasoftendifficulttoseesamplemixupsinthecsvmatrix:
hist(cornaomit,main="HistogramofRSprobes")
RS_SN<-hclust(dist(cornaomit))
#plotclustering
plot(RS_SN,labels= pheno$Sample_ID,cex=0.5,las=2,main='Hierarchial clustering 59 RS SNP probes') ###adds sample ids 



### Checking samples from unrelated individuals ~~~~~~~~~~~~~~~~~~~~~~~~~~~

brainIDs<-unique(pheno$Individual_ID) 
unrelatedcors<-snpCor
for (i in brainIDs){
  samples<-which(pheno$Individual_ID == i)
  unrelatedcors[samples,samples]<-NA
}
maxunrelatedcors<-apply(unrelatedcors, 1, max, na.rm = TRUE)
hist(maxunrelatedcors, main="Maximum correlation in samples from unrelated individuals", xlab="Max Correlation")
length(which(maxunrelatedcors>0.8))

QCmetrics<-cbind(QCmetrics,maxunrelatedcors)
SamplesFail[which(maxunrelatedcors > 0.8)]<-TRUE
Step5<-c(length(which(QCmetrics$maxunrelatedcors > 0.8)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step5)
print(Stepsummary)




### Pfilter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###The pfilter function in the wateRmelon package filters data sets based on beadcounts and detection p-values.
###Beadcounts
###If the percentage of samples with a beadcount less than 3 is greater than 5% for any probe, the probe is removed.
###Detection p-values
###If the percentage of probes with a detection p-value less than 0.05 is greater than 1% for any sample, the sample is removed.
###Similarly, if the percentage of samples with a detection p-value less than 0.05 is greater than 1% for any probe, the probe is removed.

msetEPIC.pf <- pfilter(msetEPIC)
#remove the probes that failed the pfilter
msetEPIC<-msetEPIC[rownames(betas(msetEPIC)) %in% rownames(betas(msetEPIC.pf)),]
#mark samples that fail the pfilter
pFilterPass<-colnames(betas(msetEPIC)) %in% colnames(betas(msetEPIC.pf))
QCmetrics<-cbind(QCmetrics,pFilterPass)
SamplesFail[which(pFilterPass==FALSE)]<-TRUE
Step6<-c(length(which(pFilterPass==FALSE)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step6)

### Outliers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##The outlyx function in the wateRmelon package can be used to check if any of the remaining samples are classed as 'outliers'.
betas <- betas(msetEPIC)
outliers <- outlyx(betas)
QCmetrics <- cbind(QCmetrics, outliers$outliers)
SamplesFail[which(outliers$outliers == TRUE)]<-TRUE
Step7<-c(length(which(outliers$outliers == TRUE)),sum(SamplesFail))
Stepsummary<-cbind(Stepsummary,Step7)

### Normalisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### The methylation data for just the passed samples and probes is quantile normalised using the dasen function from the wateRmelon package.
### Normalisation is done seperately for each brain region

#this removes all the failed samples from the mset (pfiltered probes already removed)
msetEPIC<-msetEPIC[,!SamplesFail]

###subset pheno by tissue type
pheno_region <-pheno_Sample_Call_Pass[pheno_Sample_Call_Pass$Tissue_Type=="region",] 
QCmetrics_region <- QCmetrics_Sample_Call[QCmetrics_Sample_Call$Tissue_Type=="region",]
dim(QCmetrics_region)
### pulls out ids from region pheno file
region_ids <- as.character(pheno_region$Basename)
length(region_ids)

###subsets the region samples from the mset and reassigns 
msetEPIC_region <- msetEPIC[,(colnames(msetEPIC)%in% region_ids)]

### run dasen normalisation
msetEPIC_region.dasen<-dasen(msetEPIC_region)

###saved the normalised betas of samples that have passed from region
save(msetEPIC_region.dasen, file="normalised_betas_region.rdat") 

plotmset_density<-function(mset, study="Density Plots HD-Study"){
  onetwo<-fData(mset)$DESIGN
  mat<-betas(mset)
  
  plot(density(mat[onetwo=="I",1], na.rm=T, bw=0.03), cex.main=0.8, main=paste(study, "Betas"), ylim=c(0, 5.2), xlab="")
  lines(density(mat[onetwo=="II",1], na.rm=T, bw=0.03), col="red")
  
  for(j in 2:ncol(mat)){
    lines(density(mat[onetwo=="I",j], na.rm=T, bw=0.03))
    lines(density(mat[onetwo=="II",j], na.rm=T, bw=0.03), col="red")
  }
  
  legend("topright", legend=c("Type I", "Type II"), lty=1, col=c("black", "red")) 
}
plotmset_density(msetEPIC_region, study="Huntington's Disease region")
plotmset_density(msetEPIC_region.dasen, study="Huntington's Disease region Dasen")


###Cell composition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###neuronal cell estimates using CETS package
load("cetsBrain.rda")
load("cetsDilution.rda")

## INDEX CELL TYPES
modelIdx <- list(neuron = pdBrain$celltype == "N", glia = pdBrain$celltype == "G")

## MAKE REF FILES
refProfile <- getReference(brain, modelIdx)
head(refProfile)
dim(refProfile)

### MATCH CPGS IN REF AND TARGET FILES
betas.dasen<-betas(msetEPIC_region.dasen)
betas.dasen2<-betas.dasen[which(rownames(betas.dasen)%in%rownames(refProfile)),]
dim(betas.dasen2) # 
refProfile2<-refProfile[which(rownames(refProfile)%in%rownames(betas.dasen2)),]
dim(refProfile2) # 

### ESTIMATE PROPORTIONS OF NEURONAL COMPOSITION 
prop <- estProportion(betas.dasen2, profile = refProfile2)
head(prop)
prop<-as.data.frame(prop)
prop<-prop[match(QCmetrics_region$Basename, rownames(prop)),]
QCmetrics_region<-cbind(QCmetrics_region, prop)


### filter probes from normalised betas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
betas<-betas(msetEPIC_region.dasen)

crosshyb<-read.table("/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
snpProbes<-read.table("SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]

betas<-betas[!(rownames(betas) %in% crosshyb[,1]), ]
betas<-betas[!(rownames(betas) %in% unique(snpProbes$IlmnID)), ]
betas<-betas[-grep("rs", rownames(betas)),]

#only keeping the useful variables in pheno
pheno_region<-QCmetrics_region
pheno_region<-subset(pheno_region, select=-c(Control,Intensity,PredictedSex,pFilterPass))

###Finally the QCed, normalised and filtered dataset is saved 
###to HD_Normalised.rdat ready for subsequent analysis.
#save QCed object
save(pheno_region, betas, file = "HD_Normalised_region.rdat")