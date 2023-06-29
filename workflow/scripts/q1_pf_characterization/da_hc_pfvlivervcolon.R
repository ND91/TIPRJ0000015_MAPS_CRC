#!/usr/bin/env Rscript

# This script will perform the differential abundance analyses comparing liver and colon with PF for the desired comparison.

if(require(speckle) == FALSE){
  devtools::install_github("phipsonlab/speckle")
  library(speckle)
} else{
  library(speckle)
}

library(dplyr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste0("Script needs 9 arguments. Current input is:", args))
}

seurat_pbmc_pf_liver_colon_rds <- args[1]
seuratDA_r <- args[2]
comparison <- args[3]
dacs_anova_rds <- args[4]
dacs_anova_csv <- args[5]
dacs_livervpf_rds <- args[6]
dacs_livervpf_csv <- args[7]
dacs_colonvpf_rds <- args[8]
dacs_colonvpf_csv <- args[9]

source(seuratDA_r)

seurat_pbmc_pf_liver_colon <- readRDS(seurat_pbmc_pf_liver_colon_rds)

seurat_pf_liver_colon_cells <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::filter(Tissue %in% c("PF", "Liver", "Colon")) %>%
  dplyr::pull(CellID)
seurat_pf_liver_colon <- seurat_pbmc_pf_liver_colon[,which(seurat_pbmc_pf_liver_colon@meta.data$CellID %in% seurat_pf_liver_colon_cells)]
sampleinfo_pf_liver_colon <- seurat_pf_liver_colon@meta.data %>%
  dplyr::select(SampleID, Donor, Tissue) %>%
  unique()

seurat_pf_liver_cells <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::filter(Tissue %in% c("PF", "Liver")) %>%
  dplyr::pull(CellID)
seurat_pf_liver <- seurat_pbmc_pf_liver_colon[,which(seurat_pbmc_pf_liver_colon@meta.data$CellID %in% seurat_pf_liver_cells)]
sampleinfo_pf_liver <- seurat_pf_liver@meta.data %>%
  dplyr::select(SampleID, Donor, Tissue) %>%
  unique()

seurat_pf_colon_cells <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::filter(Tissue %in% c("PF", "Colon")) %>%
  dplyr::pull(CellID)
seurat_pf_colon <- seurat_pbmc_pf_liver_colon[,which(seurat_pbmc_pf_liver_colon@meta.data$CellID %in% seurat_pf_colon_cells)]
sampleinfo_pf_colon <- seurat_pf_colon@meta.data %>%
  dplyr::select(SampleID, Donor, Tissue) %>%
  unique()

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

# ANOVA

if(reference_level != "manual_l0"){
  dacs_anova <- seuratDA(seuratobj = seurat_pf_liver_colon, 
                         sampleinfo = sampleinfo_pf_liver_colon, 
                         cellsampleID = "SampleID", 
                         cellclusterID = case_level, 
                         clusterparentID = reference_level, 
                         ncell_threshold = 20, 
                         method = "propeller", 
                         design = "Tissue", 
                         return_props = T)
  
  dacs_anova_df <- do.call(rbind, lapply(names(dacs_anova), function(celltype){
    celltype_df <- dacs_anova[[celltype]]$da
    celltype_df$reference <- celltype
    return(celltype_df)
  }))
  
  dacs_livervpf <- seuratDA(seuratobj = seurat_pf_liver, 
                            sampleinfo = sampleinfo_pf_liver, 
                            cellsampleID = "SampleID", 
                            cellclusterID = case_level, 
                            clusterparentID = reference_level, 
                            ncell_threshold = 20, 
                            method = "propeller", 
                            design = "Tissue", 
                            return_props = T)
  
  dacs_livervpf_df <- do.call(rbind, lapply(names(dacs_livervpf), function(celltype){
    celltype_df <- dacs_livervpf[[celltype]]$da
    celltype_df$reference <- celltype
    return(celltype_df)
  }))
  
  dacs_colonvpf <- seuratDA(seuratobj = seurat_pf_colon, 
                            sampleinfo = sampleinfo_pf_colon, 
                            cellsampleID = "SampleID", 
                            cellclusterID = case_level, 
                            clusterparentID = reference_level, 
                            ncell_threshold = 20, 
                            method = "propeller", 
                            design = "Tissue", 
                            return_props = T)
  
  dacs_colonvpf_df <- do.call(rbind, lapply(names(dacs_colonvpf), function(celltype){
    celltype_df <- dacs_colonvpf[[celltype]]$da
    celltype_df$reference <- celltype
    return(celltype_df)
  }))
  
} else{
  dacs_anova <- seuratDA(seuratobj = seurat_pf_liver_colon, 
                         sampleinfo = sampleinfo, 
                         cellsampleID = "SampleID", 
                         cellclusterID = case_level, 
                         ncell_threshold = 20, 
                         method = "propeller", 
                         design = "Tissue",
                         return_props = T)
  dacs_anova[[1]]$da$reference <- "All"
  dacs_anova_df <- dacs_anova[[1]]$da
  
  dacs_livervpf <- seuratDA(seuratobj = seurat_pf_liver, 
                            sampleinfo = sampleinfo_pf_liver, 
                            cellsampleID = "SampleID", 
                            cellclusterID = case_level, 
                            ncell_threshold = 20, 
                            method = "propeller", 
                            design = "Tissue",
                            return_props = T)
  dacs_livervpf[[1]]$da$reference <- "All"
  dacs_livervpf_df <- dacs_livervpf[[1]]$da
  
  dacs_colonvpf <- seuratDA(seuratobj = seurat_pf_colon, 
                            sampleinfo = sampleinfo_pf_colon, 
                            cellsampleID = "SampleID", 
                            cellclusterID = case_level, 
                            ncell_threshold = 20, 
                            method = "propeller", 
                            design = "Tissue",
                            return_props = T)
  dacs_colonvpf[[1]]$da$reference <- "All"
  dacs_colonvpf_df <- dacs_colonvpf[[1]]$da
}

saveRDS(dacs_anova, dacs_anova_rds, compress = "gzip")
write.csv(dacs_anova_df, dacs_anova_csv)
saveRDS(dacs_livervpf, dacs_livervpf_rds, compress = "gzip")
write.csv(dacs_livervpf_df, dacs_livervpf_csv)
saveRDS(dacs_colonvpf, dacs_colonvpf_rds, compress = "gzip")
write.csv(dacs_colonvpf_df, dacs_colonvpf_csv)

sessionInfo()