### load libraries 

library("tidyverse") 
library("EmpiricalBrownsMethod")


load("model.RData") # load annotated EWAS model results 
load("Normalised_betas.rdat") # load normalised betas


### Create a matrix of genomic regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### associated with HD AOO according to Lee et al 2019including 35kb upstream and 10kb downstream

GRegions <- matrix(c("chr15",    "31185854", "31293807", "MTMR10",
"chr15",    "31156703", "31245311", "FAN1",
"chr15",    "31258264", "31463476", "TRPM1 ",
"chr4", "2810454",  "2941803",  "ADD1",
"chr4", "2592159",  "2744302",  "FAM193A",
"chr4", "2904663",  "2975233",  "NOP14",
"chr4", "2236324",  "2430370",  "ZFYVE28",
"chr3", "36992357", "37044795", "EPM2AIP1",
"chr19",    "48516100", "48624109", "PLA2G4C",
"chr4", "2759750",  "2852823",  "SH3BP2",
"chr3", "37059117", "37227992", "LRRFIP2",
"chr19",    "48583702", "48683852", "LIG1",
"chr4", "1778206",  "1867974",  "LETM1",
"chr3", "36999841", "37102337", "MLH1",
"chr8", "103181729",    "103261346",    "RRM2B",
"chr5", "79915467", "80182634", "MSH3",
"chr4", "1659527",  "1724421",  "SLBP",
"chr4", "1838123",  "1993934",  "WHSC1",
"chr4", "2038645",  "2253848",  "POLN",
"chr19",    "48638952", "48710879", "C19orf68",
"chr4", "1688217",  "1756905",  "TACC3",
"chr5", "79887045", "79960800", "DHFR",
"chr7", "5930777",  "6020322",  "RSPH10B",
"chr4", "2930232",  "3052474",  "GRK4",
"chr5", "79817574", "79876304", "ANKRD34B",
"chr4", "1606608",  "1695988",  "FAM53A",
"chr2", "190504018",    "190635919",    "ANKAR",
"chr2", "190576386",    "190638020",    "OSGEPL1",
"chr7", "5977870",  "6058737",  "PMS2",
"chr5", "79910819", "79956854", "MTRNR2L2",
"chr4", "2897288",  "2946586",  "MFSD10",
"chr7", "6013882",  "6073465",  "AIMP2",
"chr7", "6026878",  "6108860",  "EIF2AK1",
"chr4", "3041408",  "3255687",  "HTT",
"chr4", "2708387",  "2768103",  "TNIP2",
"chr2", "190599993",    "190659097",    "ORMDL1 ",
"chr7", "6036007",  "6086183",  "ANKRD61",
"chr3", "37391477", "37486988", "C3orf35",
"chr15",    "30865445", "30916677", "LOC101927579",
"chr3", "37249682", "37418370", "GOLGA4",
"chr15",    "30883879", "30941013", "ARHGAP11B",
"chr2", "190491125",    "190545557",    "ASNSD1",
"chr2", "190613710",    "190752355",    "PMS1",
"chr5", "145791873",    "145901071",    "TCERG1  ",
"chr5", "145859417",    "145905676",    "GPR151"
),ncol = 4, byrow = T)

### Create BEDfile of genomic regions of interest 

write.table(GRegions, file="GRegions.bed", sep="\t", col.names= FALSE,row.names=F, quote=FALSE )


### create BEDfile of EPIC probes  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GR.Probes<-data.frame(paste0("chr",model$CHR), 
                   as.numeric(model$MAPINFO),
                   as.numeric(model$MAPINFO),
                   modelCer$IlmnID,
                   as.numeric(modelCer$pval.bc))

GR.Probes <- na.omit(GR.Probes)

write.table(GR.Probes, file="GR.probes.Cer.bed", sep="\t", col.names= FALSE,row.names=F, quote=FALSE )



### Identification of the EPIC probes falling into the defined LD regions run  on the command line ~~~~~~~~~~~~~~~~~~

# ml BEDTools/2.27.1-foss-2018b  ### load bedtools 

# bedtools sort -i GRegions.bed > GRegions.sorted.bed
# bedtools sort -i GR.probes.bed > GR.probes.sorted.bed

### intersect the EPIC probes within the genomic regions 
# bedtools intersect -a GR.probes.sorted.bed -b GRegions.sorted.bed -wo > GR.bed

GR <- read.table("GR.bed",sep="\t",stringsAsFactors=FALSE)



#### pull out the individual genomoic regions with multiple probes on the EPIC array

regs <- GR %>% select(V9) %>% group_by(V9) %>% mutate(count = n()) %>% unique()
regs <- as.data.frame(regs)  
regs <- regs$V9

#x <- "FAN1"
regs <- as.data.frame(regs)
regs[,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- NA

regs$regs-> rownames(regs)

### for loop to run empiricial Brown's method on each region ### 
for(x in rownames(regs)){

    probes <- GRCer[GRCer$V9 %in% x, ] ### pull out probes for each gene

    N <-nrow(probes)#### get the number  
    CHR <- probes[1,1]  
    start <- probes[1,7]
    end <- probes [1,8]
    ### create pval vector of the probes 
    p.val <- probes$V5

    ### pull out betas of the probes 
    probes_betas <- betas[rownames(betas) %in% probes$V4, ]
    probes_betas<- probes_betas[match(probes$V4,rownames(probes_betas)),] ### get in the same order 


    ### run browns

    res <- empiricalBrownsMethod(probes_betas, p.val, extra_info=T)

    c(unlist(res),N, CHR, start, end) -> regs[x, c("P_test","P_Fisher", "Scale_Factor_C","DF", "N", "CHR", "Locus.start", "Locus.end")]
 
}


write.csv(regs, file = "browns_enrichment.csv")

