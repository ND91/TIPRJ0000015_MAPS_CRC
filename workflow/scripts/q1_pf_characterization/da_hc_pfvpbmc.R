#!/usr/bin/env Rscript

# This script will perform the differential abundance analyses comparing PF (case) with PBMC (reference) for the desired comparison.

if(require(speckle) == FALSE){
  devtools::install_github("phipsonlab/speckle")
  library(speckle)
} else{
  library(speckle)
}

library(dplyr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seuratDA_r <- args[2]
comparison <- args[3]
dacs_rds <- args[4]
dacs_csv <- args[5]

source(seuratDA_r)

seuratObject <- readRDS(seurat_rds)

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

sampleinfo <- seuratObject@meta.data %>%
  dplyr::select(SampleID, Donor, Tissue) %>%
  unique()

if(reference_level != "manual_l0"){
  dacs <- seuratDA(seuratobj = seuratObject, 
                   sampleinfo = sampleinfo, 
                   cellsampleID = "SampleID", 
                   cellclusterID = case_level, 
                   clusterparentID = reference_level, 
                   ncell_threshold = 20, 
                   method = "propeller", 
                   design = "Tissue", 
                   return_props = T)
  
  dacs_df <- do.call(rbind, lapply(names(dacs), function(celltype){
    celltype_df <- dacs[[celltype]]$da
    celltype_df$reference <- celltype
    return(celltype_df)
  }))
  
} else{
  dacs <- seuratDA(seuratobj = seuratObject, 
                   sampleinfo = sampleinfo, 
                   cellsampleID = "SampleID", 
                   cellclusterID = case_level, 
                   ncell_threshold = 20, 
                   method = "propeller", 
                   design = "Tissue",
                   return_props = T)
  
  dacs[[1]]$da$reference <- "All"
  
  dacs_df <- dacs[[1]]$da
}

saveRDS(dacs, dacs_rds, compress = "gzip")
write.csv(dacs_df, dacs_csv)

sessionInfo()