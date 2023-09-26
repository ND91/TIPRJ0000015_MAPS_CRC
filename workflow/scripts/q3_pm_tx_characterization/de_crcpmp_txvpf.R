#!/usr/bin/env Rscript

# This script will perform the differential expression analyses comparing PF (case) with PBMC (reference) for cells present in both only.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

crcpmp_pf_tx_seuratobject_rds <- args[1]
seuratDE_r <- args[2]
celltype_level <- args[3]
deseq2_list_rds <- args[4]
txvpf_degs_xlsx <- args[5]

source(seuratDE_r)

seuratObject <- readRDS(crcpmp_pf_tx_seuratobject_rds)

tissue_celltype <- unique(seuratObject@meta.data[,c("Tissue", celltype_level)])
celltypes_present <- names(which(colSums(table(tissue_celltype$Tissue, tissue_celltype[,celltype_level]) != 0*1) == 2)) #Identify which celltypes are present in both tissues
seuratObject_filtered <- seuratObject[,seuratObject@meta.data[,celltype_level] %in% celltypes_present]

pb_sample_metadata <- unique(seuratObject_filtered@meta.data[,c("SampleID", "Donor", "Tissue", "Sex", "Age")])
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

degs_txvpf <- seuratDE(seuratobj = seuratObject_filtered, 
                       cellsampleID = "SampleID", 
                       cellclusterID = celltype_level, 
                       sampleinfo = pb_sample_metadata, 
                       design = "~Tissue+Donor", 
                       contrast = c("Tissue", "TX", "PF"))

degs_txvpf <- degs_txvpf[lapply(degs_txvpf, length) != 0]

saveRDS(degs_txvpf, deseq2_list_rds, compress = "gzip")

names(degs_txvpf) <- substr(names(degs_txvpf), 1, 30)

write.xlsx(lapply(degs_txvpf, function(celltype){
  celltype$degs
}), file = txvpf_degs_xlsx)

sessionInfo()
