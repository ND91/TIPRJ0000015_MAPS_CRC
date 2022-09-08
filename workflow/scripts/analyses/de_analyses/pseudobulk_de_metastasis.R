#!/usr/bin/env Rscript
# This script will perform pseudobulk DE analysis comparing differences between metastatic and non-metastatic samples using DESeq2.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
deseq2_list_path <- args[2]
celltype_level <- args[3]
tissue <- args[4]
functions_r_path <- args[5]

source(functions_r_path)

seuratObject <- readRDS(seurat_rds_path)

seuratObject <- seuratObject[,seuratObject@meta.data$Tissue == tissue]

sampleinfo <- unique(seuratObject@meta.data[,c("Donor", "Sex", "Age", "Group")])
rownames(sampleinfo) <- paste0(sampleinfo$Donor)
sampleinfo$Metastasis <- ifelse(sampleinfo$Group == "CRC+", "PMp", "PMn")

deseq2_results_dds_list <- seuratDE(seuratobj = seuratObject, 
                                    cellsampleID = "Donor", 
                                    cellclusterID = celltype_level, 
                                    sampleinfo = sampleinfo, 
                                    design = "~Metastasis+Sex+Age", 
                                    contrast = c("Metastasis", "PMp", "PMn"))

saveRDS(deseq2_results_dds_list, deseq2_list_path, compress = "gzip")

sessionInfo()