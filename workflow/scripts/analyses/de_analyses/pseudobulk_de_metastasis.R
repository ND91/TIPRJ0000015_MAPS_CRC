#!/usr/bin/env Rscript
# This script will perform pseudobulk DE analysis comparing differences between metastatic and non-metastatic samples using DESeq2.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
degs_list_path <- args[2]
deseq2_list_path <- args[3]
base_path <- args[4]
celltype_level <- args[5]
tissue <- args[6]
functions_r_path <- args[7]

source(functions_r_path)

seuratObject <- readRDS(seurat_rds_path)

seuratObject <- seuratObject[,seuratObject@meta.data$Tissue == tissue]

sampleinfo <- unique(seuratObject@meta.data[,c("Donor", "Sex", "Age", "Metastasis")])
rownames(sampleinfo) <- paste0(sampleinfo$Donor)
sampleinfo$Metastasis <- ifelse(sampleinfo$Metastasis == "PM+", "PMp", "PMn")

deseq2_results_dds_list <- seuratDE(seuratobj = seuratObject, 
                                    cellsampleID = "Donor", 
                                    cellclusterID = celltype_level, 
                                    sampleinfo = sampleinfo, 
                                    design = "~Metastasis", 
                                    contrast = c("Metastasis", "PMp", "PMn"))

dds_list <- degs_list <- vector(mode = "list", length = length(deseq2_results_dds_list))
names(dds_list) <- names(degs_list) <- names(deseq2_results_dds_list)

for(celltype in names(deseq2_results_dds_list)){
  degs_list[[celltype]] <- deseq2_results_dds_list[[celltype]]$degs
  dds_list[[celltype]] <- deseq2_results_dds_list[[celltype]]$dds
  
  celltype_rn <- gsub("\\+", "p", celltype)
  celltype_rn <- gsub("[\\| \\/]", "_", celltype_rn)
  
  write.csv(degs_list[[celltype]], paste0(base_path, "_", celltype_rn, "_degs.csv"))
  saveRDS(dds_list[[celltype]], paste0(base_path, "_", celltype_rn, "_dds.Rds"))
}

saveRDS(degs_list, degs_list_path, compress = "gzip")
saveRDS(dds_list, deseq2_list_path, compress = "gzip")

sessionInfo()