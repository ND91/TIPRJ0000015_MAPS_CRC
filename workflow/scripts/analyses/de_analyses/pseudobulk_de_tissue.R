#!/usr/bin/env Rscript
# This script will perform pseudobulk DE analysis comparing diferences between PF, PBMC, and TX using DESeq2.

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
group_tissue <- unlist(strsplit(args[6], "_")) 
tissue_comparison <- unlist(strsplit(group_tissue[2], "v")) #This is why we need to use the "v" as separator
functions_r_path <- args[7]

source(functions_r_path)

tissue_case <- tissue_comparison[1]
tissue_control <- tissue_comparison[2]

seuratObject <- readRDS(seurat_rds_path)

sampleinfo <- unique(seuratObject@meta.data[,c("Donor", "Tissue", "Sex", "Age")])
rownames(sampleinfo) <- paste0(sampleinfo$Donor, "_", sampleinfo$Tissue)

group <- ifelse(group_tissue[1] == "PMp", "PM+", "PM-")

seuratObject <- seuratObject[,seuratObject@meta.data$Metastasis == group]

deseq2_results_dds_list <- seuratDE(seuratobj = seuratObject, 
                                    cellsampleID = "SampleID", 
                                    cellclusterID = celltype_level, 
                                    sampleinfo = sampleinfo, 
                                    design = "~Tissue", 
                                    contrast = c("Tissue", tissue_case, tissue_control))

degs_list <- lapply(names(deseq2_results_dds_list), function(celltype){
  degs <- deseq2_results_dds_list[[celltype]]$degs

  celltype_rn <- gsub("\\+", "p", celltype)
  celltype_rn <- gsub("[\\| \\/]", "_", celltype_rn)
  
  write.csv(degs, paste0(base_path, "_", celltype_rn, ".csv"))
  
  return(degs)
})

saveRDS(degs_list, degs_list_path, compress = "gzip")

dds_list <- lapply(deseq2_results_dds_list, function(celltype){
  celltype$dds
})
names(dds_list) <- names(deseq2_results_dds_list)

saveRDS(dds_list, deseq2_list_path, compress = "gzip")

sessionInfo()