#!/usr/bin/env Rscript
# This script will perform DA analysis comparing differences between metastatic with non-metastatic using propeller.

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
tissue <- args[4]
functions_r_path <- args[5]

source(functions_r_path)

level_child <- relative[1]
level_parent <- relative[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject@meta.data$Metastasis <- ifelse(seuratObject@meta.data$Metastasis == "PM+", "PMp", "PMn")

seuratObject <- seuratObject[,seuratObject@meta.data$Tissue == tissue]

if(level_parent == "NULL"){
  propeller_list <- seuratDA(seuratobj = seuratObject, 
                             cellsampleID = "SampleID", 
                             cellclusterID = level_child,
                             method = "propeller", 
                             design = "Metastasis")  
} else{
  propeller_list <- seuratDA(seuratobj = seuratObject, 
                             cellsampleID = "SampleID", 
                             clusterparentID = level_parent,
                             cellclusterID = level_child,
                             method = "propeller", 
                             design = "Metastasis")  
}

saveRDS(propeller_list, propeller_list_path, compress = "gzip")