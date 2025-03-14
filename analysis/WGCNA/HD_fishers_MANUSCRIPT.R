
library("tidyverse") 

### load file with probes in module of interest 
col_probes <- read.csv("MMcol.csv", stringsAsFactors = FALSE, row.names = 1)

### load regions associated with HD AOO according to Lee et al 2019 including 35kb upstream and 10kb downstream
HD_regions <-read.table(file="GR.bed", ,sep="\t",stringsAsFactors=FALSE)

col_probes <- col_probes$Name

regs <- HD_regions %>% select(V9) %>% group_by(V9) %>% mutate(count = n()) %>% unique()
regs <- as.data.frame(regs)
regs <- regs$V9

regs <- as.data.frame(regs)
regs[,c("region_probes","shared_probes", "p_value","shared_CpGs", "odds_ratio")] <- NA
regs$regs-> rownames(regs)

N <- 483458 #### background (nummber of probes in WGCNA dataset)
listA <- col_probes

#x <- "FAN1"

### for loop to run fisher test on each region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(x in rownames(regs)){

    probes <- HD_regions[HD_regions$V9 %in% x, ] ### pull out probes for each gene

listB <- probes$V4
  listA <- unique(listA[!is.na(listA)])
  listB <- unique(listB[!is.na(listB)])

   n <- length(listA[!is.na(listA)])
      m <- length(listB[!is.na(listB)])
      shared.CpGs <- intersect(listA , listB)
      shared.CpGs <- shared.CpGs[!is.na(shared.CpGs)]
      k <- length(shared.CpGs)
      
      con.table <- matrix(c(k, m-k, n-k, N-n-m+k),
                          nrow = 2,
                          dimnames = list(B=c("inB","notB"),
                                          A=c("inA","notA")))
      f.test <- fisher.test(con.table, alternative = "greater")
      results <- list(N.listB=m , N.Shared=k , P.value = f.test$p.value , Shared.CpGs = paste(shared.CpGs,collapse = ";"), estimate = f.test$estimate)

      unlist(results) -> regs[x, c("region_probes","shared_probes", "p_value","shared_CpGs", "odds_ratio")]

}

write.csv(regs, file = "col_HD_enrichment.csv")
