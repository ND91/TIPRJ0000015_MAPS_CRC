#!/usr/bin/env Rscript
# This script will perform pseudobulk DE analysis comparing diferences between PF, PBMC, and TX using DESeq2.

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
tissue_comparison <- unlist(strsplit(args[4], "v")) #This is why we need to use the "v" as separator
functions_r_path <- args[5]

source(functions_r_path)

tissue_case <- tissue_comparison[1]
tissue_control <- tissue_comparison[2]

seuratObject <- readRDS(seurat_rds_path)

sampleinfo <- unique(seuratObject@meta.data[,c("Donor", "Tissue", "Sex", "Age")])
rownames(sampleinfo) <- paste0(sampleinfo$Donor, "_", sampleinfo$Tissue)

deseq2_results_dds_list <- seuratDE(seuratobj = seuratObject, cellsampleID = "SampleID", cellclusterID = celltype_level, sampleinfo = sampleinfo, design = "~Tissue+Sex+Age", contrast = c("Tissue", tissue_case, tissue_control))

saveRDS(deseq2_results_dds_list, deseq2_list_path, compress = "gzip")

sessionInfo()