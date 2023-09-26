#!/usr/bin/env Rscript

# This script will perform the differential expression analyses comparing CRC with HC for the desired tissue.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seuratDE_r <- args[2]
celltype_level <- args[3]
degs_list_rds <- args[4]
degs_xlsx <- args[5]

source(seuratDE_r)

seuratObject <- readRDS(seurat_rds)

pb_sample_metadata <- unique(seuratObject@meta.data[,c("SampleID", "Donor", "Group", "Sex", "Age", "PCI")])
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

degs_pci <- seuratDE(seuratobj = seuratObject, 
                     cellsampleID = "SampleID", 
                     cellclusterID = celltype_level, 
                     sampleinfo = pb_sample_metadata, 
                     design = "~PCI",
                     name = "PCI")

degs_pci <- degs_pci[lapply(degs_pci, length) != 0]

saveRDS(degs_pci, degs_list_rds, compress = "gzip")

write.xlsx(lapply(degs_pci, function(celltype){
  celltype$degs
}), file = degs_xlsx)

sessionInfo()
