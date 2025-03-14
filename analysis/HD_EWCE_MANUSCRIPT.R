### load libraries 

library(DropletUtils)
library(EWCE)
library(ewceData)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(dplyr) 

### load barcodes from selected single cell dataset - Lee et al 2020  https://pubmed.ncbi.nlm.nih.gov/32681824/
barcodes <- read.table("file_path/GSE152058_data/barcodes.tsv.gz")
barcodes <- barcodes[-1, ] ### removed the headers from barcode 
write.table(barcodes, file ="barcodes.tsv.gz", col.names = FALSE, row.names = FALSE, quote = FALSE ) #### resave (after deleting orginal barcode file) 
																								#ensuring no row or col names are forced and no charachter vectors are force

#### repeat with features
features <- read.table("file_path/GSE152058_data/features.tsv.gz")
### removed the colnames from features as it is resaved 

write.table(features , file ="file_path/GSE152058_data/features.tsv.gz", col.names = FALSE, row.names = TRUE, quote = FALSE ) #### resave (after deleting orginal barcode file) 
																								#ensuring no col names are forced and no charachter vectors are force (keep rownames)

data_dir <- paste("/GSE152058_data")
list.files(data_dir)

### Seurat loads in the sparse matrices
sce <- Seurat::Read10X(data_dir)

### nead the cell annotation in seperately to assign colData
umapMeta <- read.table("file_path/GSE152058_human_snRNA_processed_coldata.tsv", header = T) 


### create summarized experiment object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SE <- SummarizedExperiment(assays=list(counts=sce), colData = umapMeta)

save(SE, file = "HD_cellEnrich_SE_object.Rdata")

### remove uniformative genes - for mouse data change input species to "mouse"
Expdrop <- drop_uninformative_genes(exp=SE, dge_method = "limma", drop_nonhuman_genes = T, 
										input_species = "human", level2annot = SE$CellType, no_cores = 20)


annotLevels = list(level1class = SE$CellType,
					level2class = SE$CellType)

fNames_ <- generate_celltype_data(exp = Expdrop,
										annotLevels = annotLevels,
										groupName = "sc",
										savePath = "file_path", no_cores = 20)

load("file_path/ctd_sc.rda")

### create CTD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ctd <- EWCE::load_rdata(fNames)


#### check marker gene expression with plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try({
  pltexp <- EWCE::plot_ctd(ctd = ctd,
                        level = 1,
                        genes = c("RBFOX3", "FOXP2","PPP1R1B", "DRD1","DRD2", "CHAT","GFAP","CSF1R", "OLIG2"),
                        metric = "expression")
})


try({
  pltspec <- EWCE::plot_ctd(ctd = ctd,
                        level = 1,
                        genes = c("RBFOX3", "FOXP2","PPP1R1B", "DRD1","DRD2", "CHAT","GFAP","CSF1R", "OLIG2"),
                        metric = "specificity")
})


#### read in the genelist

col_probes <- read.table(file = "col_probes.txt")
genelist <-col_probes$UCSC1

### run bootstrap enrichment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#for mouse data change sctspecies and genelistSpecies to "mouse"

unconditional_results_col_probes <- EWCE::bootstrap_enrichment_test(
	sct_data = ctd,
	hits = genelist,
	sctSpecies = "human",
	genelistSpecies = "human",
	reps = 100000,
	annotLevel = 1,
	no_cores = 20)


### plot enrichment

cellRes <- unconditional_col_probes$results

cell_order <- c("T_Cell", "Mural", "Endothelial","PV_Interneuron","GABAergic_Interneuron","Cholinergic_Interneuron",
				 "FOXP2_Neur","D2_MSN", "D1_MSN", "Sec_Ependymal", "Cil_Ependymal","OPC","Oligodendrocyte","Microglia","Astrocyte")

cellRes <- cellRes %>%
  slice(match(cell_order, CellType))

cellRes$CellType <- factor(cellRes$CellType, levels = cell_order)

# Remove negative sd_from_mean for plotting (as in the EWCE plotting function) 
cellRes[which(cellRes$sd_from_mean < 0),"sd_from_mean"] <- NA

#Define significant enrichments 
cellRes[which(cellRes$q < 0.05),"sig"] <- "  *"
cellRes[which(cellRes$q > 0.05),"sig"] <- ""

  ggplot(cellRes, aes(x = as.factor(CellType), y = sd_from_mean))+
  geom_bar(stat = "identity", color = "black",fill = "red")+
  geom_text(aes(label = sig, x = as.factor(CellType), y=sd_from_mean + .05), size = 10)+
  xlab(NULL)+
  ylab("Standard Deviations from Mean")+
  coord_flip()+
  theme_bw(base_size = 22)+ 
             theme(legend.position = "none", 
                   plot.margin = margin(r = 0.1,
                                        l = 0.1),
                   plot.title = element_text(hjust = 0.5, vjust = 0.5))+
             ylim(0, 13)

