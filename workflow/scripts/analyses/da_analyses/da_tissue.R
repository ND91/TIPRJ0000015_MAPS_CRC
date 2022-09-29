#!/usr/bin/env Rscript
# This script will perform DA analysis comparing differences between PF, PBMC, and TX using propeller.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
propeller_list_path <- args[2]
relative <- unlist(strsplit(args[3], "r")) #This is why we need to use the "r" as separator
group_tissue <- unlist(strsplit(args[4], "_")) 
tissue_comparison <- unlist(strsplit(group_tissue[2], "v")) #This is why we need to use the "v" as separator
functions_r_path <- args[5]

source(functions_r_path)

level_child <- relative[1]
level_parent <- relative[2]

tissue_case <- tissue_comparison[1]
tissue_control <- tissue_comparison[2]

seuratObject <- readRDS(seurat_rds_path)

group <- ifelse(group_tissue[1] == "PMp", "PM+", "PM-")

seuratObject <- seuratObject[,seuratObject@meta.data$Tissue %in% c(tissue_case, tissue_control) & seuratObject@meta.data$Metastasis == group]

if(level_parent == "NULL"){
  propeller_list <- seuratDA(seuratobj = seuratObject, 
                             cellsampleID = "SampleID", 
                             cellclusterID = level_child,
                             method = "propeller", 
                             design = "Tissue")  
} else{
  propeller_list <- seuratDA(seuratobj = seuratObject, 
                             cellsampleID = "SampleID", 
                             clusterparentID = level_parent,
                             cellclusterID = level_child,
                             method = "propeller", 
                             design = "Tissue")  
}

saveRDS(propeller_list, propeller_list_path, compress = "gzip")