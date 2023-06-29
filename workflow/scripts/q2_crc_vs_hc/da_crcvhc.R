#!/usr/bin/env Rscript

# This script will perform the differential abundance analyses comparing CRC with HC for the desired tissue.

if(require(speckle) == FALSE){
  devtools::install_github("phipsonlab/speckle")
  library(speckle)
} else{
  library(speckle)
}

library(dplyr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seuratDA_r <- args[2]
tissue <- args[3]
comparison <- args[4]
dacs_rds <- args[5]
dacs_csv <- args[6]
proportions_csv <- args[7]

source(seuratDA_r)

seuratObject <- readRDS(seurat_rds)

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

sampleinfo <- seuratObject@meta.data %>%
  dplyr::select(SampleID, Donor, Tissue, Group) %>%
  dplyr::filter(Tissue == tissue) %>%
  unique()

if(reference_level != "manual_l0"){
  dacs <- seuratDA(seuratobj = seuratObject, 
                   sampleinfo = sampleinfo, 
                   cellsampleID = "SampleID", 
                   cellclusterID = case_level, 
                   clusterparentID = reference_level, 
                   ncell_threshold = 20, 
                   method = "propeller", 
                   design = "Group", 
                   return_props = T)
  
  dacs_filtered <- dacs[which(unlist(lapply(dacs, function(celltype){length(celltype$da)})) != 0)]
  
  dacs_df <- do.call(rbind, lapply(names(dacs_filtered), function(celltype){
    celltype_df <- dacs_filtered[[celltype]]$da
    celltype_df$reference <- celltype
    return(celltype_df)
  }))
  
  proportions_df <- do.call(rbind, lapply(names(dacs_filtered), function(celltype){
    celltype_df <- data.frame(dacs_filtered[[celltype]]$props$Proportions)
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
                   design = "Group",
                   return_props = T)
  
  dacs[[1]]$da$reference <- "All"
  
  dacs_df <- dacs[[1]]$da
  
  proportions_df <- data.frame(dacs[[1]]$props$Proportions)
}

saveRDS(dacs, dacs_rds, compress = "gzip")
write.csv(dacs_df, dacs_csv)
write.csv(proportions_df, proportions_csv)

sessionInfo()