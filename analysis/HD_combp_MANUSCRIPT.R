
### load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(wateRmelon)
library(methylumi)
library(dplyr)
library(qqman)


load("model.RData")

######## Create files for comb-p ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###creates function
DMRinput <- function(results,  chr = "CHR", 
                     mapinfo = "MAPINFO", p = "pval.bc"){
  dmr <- data.frame("chrom" = paste0("chr", as.character(results[ , chr])),
                    "start" = results[, mapinfo],
                    "end" = results[, mapinfo] + 1,
                    "pvalue" = results[ , p])
  dmr <- summarise(group_by(dmr, chrom, start), end = end, pvalue = pvalue)
  colnames(dmr) <- c("#chrom", "start", "end", "pvalue")
  return(dmr)
}

###creates dataframe
CvHD<-DMRinput(model, p="pval.bc")

CvHD<-na.omit(CvHD)
summary(CvHD)

###writes as tabel for input into pipeline 
write.table(CvHD, file = "ControlvHD_DMR.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)

#run in comb-p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### this step is run through the commandline 
#   comb-p pipeline -c 4 --dist 500 --seed 1e-2 --anno hg19 -p /filepath/ControlvHD_DMR /filepath/ControlvHD_DMR.bed

### return to R

col<-unlist(strsplit(readLines("/filepath/ControlvHD_DMR.anno.hg19.bed", n = 1), split = "\t"))
col[1]<-"chr"
####adds names of columns 

###reads in annotated dmrs
dmrs <- read.table("/filepath/ControlvHD_DMR.anno.hg19.bed",
                      sep="\t",stringsAsFactors=FALSE,col.names=col)

#### sidak corrected DMRs
dmrs <- dmrs[dmrs$z_sidak_p < 0.05 & dmrs$n_probes>2,]
dmrs<-dmrs[order(dmrs$z_sidak_p),]


model[ , c("CHR","MAPINFO")]->subset 
probesdmrs <- data.frame()
for(i in 1:nrow(dmrs)){
  probesdmrs <- rbind(probesdmrs, cbind(subset[subset$CHR== as.integer(substring(dmrs$chr[i],4)) &subset$MAPINFO >= dmrs$start[i] & subset$MAPINFO <= dmrs$end[i], ]))
}

probesdmrs <- na.omit(probesdmrs)


### plot with miniman function from https://github.com/UoE-Dementia-Genomics/Created_Functions/blob/main/miniman.r
source("/filepath/miniman.R")

miniman(data = model, chr = chr, range = xxx:xxx, result = "pval.bc", 
        pad = 30000, ESlevel = 0.01, ESdat = "es.b", negcol = "red", 
        poscol = "forestgreen", cexgene = 0.6, multiply = NULL, cpgcol ="forestgreen",
        sigline = FALSE, siglinecol = "red", siglevel = siglevel,chrcol = "CHR")
