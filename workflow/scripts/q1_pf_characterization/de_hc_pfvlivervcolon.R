#!/usr/bin/env Rscript

# This script will perform the differential expression analyses comparing PF (case) with PBMC (reference) for cells present in both only.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

hc_pbmc_pf_seuratobject_rds <- args[1]
seuratDE_r <- args[2]
deseq2_pfvliver_list_rds <- args[3]
degs_pfvliver_xlsx <- args[4]
deseq2_pfvcolon_list_rds <- args[5]
degs_pfvcolon_xlsx <- args[6]

source(seuratDE_r)

seuratObject <- readRDS(hc_pbmc_pf_seuratobject_rds)

# PFvLiver
pf_liver_seuratObject <- seuratObject[,seuratObject@meta.data$Tissue %in% c("PF", "Liver")]
pfvliver_celltype <- unique(pf_liver_seuratObject@meta.data[,c("Tissue", "manual_l3")])
pfvliver_celltypes_present <- names(which(colSums(table(pfvliver_celltype$Tissue, pfvliver_celltype[,"manual_l3"]) != 0*1) >= 2)) #Identify which celltypes are present in both tissues
pf_liver_filtered_seuratObject <- pf_liver_seuratObject[,pf_liver_seuratObject@meta.data[,"manual_l3"] %in% pfvliver_celltypes_present]

pb_pfvliver_sample_metadata <- unique(pf_liver_filtered_seuratObject@meta.data[,c("SampleID", "Donor", "Tissue", "Sex", "Age")])
rownames(pb_pfvliver_sample_metadata) <- pb_pfvliver_sample_metadata$SampleID

pfvliver_degs <- seuratDE(seuratobj = pf_liver_filtered_seuratObject, 
                          cellsampleID = "SampleID", 
                          cellclusterID = "manual_l3", 
                          sampleinfo = pb_pfvliver_sample_metadata, 
                          design = "~Tissue", 
                          contrast = c("Tissue", "PF", "Liver"))

saveRDS(pfvliver_degs, deseq2_list_rds, compress = "gzip")

write.xlsx(lapply(pfvliver_degs, function(celltype){
  celltype$degs
}), file = degs_pfvliver_xlsx)

# PFvColon

pf_colon_seuratObject <- seuratObject[,seuratObject@meta.data$Tissue %in% c("PF", "Colon")]
pfvcolon_celltype <- unique(pf_colon_seuratObject@meta.data[,c("Tissue", "manual_l3")])
pfvcolon_celltypes_present <- names(which(colSums(table(pfvcolon_celltype$Tissue, pfvcolon_celltype[,"manual_l3"]) != 0*1) >= 2)) #Identify which celltypes are present in both tissues
pf_colon_filtered_seuratObject <- pf_colon_seuratObject[,pf_colon_seuratObject@meta.data[,"manual_l3"] %in% pfvcolon_celltypes_present]

pb_pfvcolon_sample_metadata <- unique(pf_colon_filtered_seuratObject@meta.data[,c("SampleID", "Donor", "Tissue", "Sex", "Age")])
rownames(pb_pfvcolon_sample_metadata) <- pb_pfvcolon_sample_metadata$SampleID

pfvcolon_degs <- seuratDE(seuratobj = pf_colon_filtered_seuratObject, 
                          cellsampleID = "SampleID", 
                          cellclusterID = "manual_l3", 
                          sampleinfo = pb_pfvcolon_sample_metadata, 
                          design = "~Tissue", 
                          contrast = c("Tissue", "PF", "Colon"))

saveRDS(pfvcolon_degs, deseq2_list_rds, compress = "gzip")

write.xlsx(lapply(pfvcolon_degs, function(celltype){
  celltype$degs
}), file = degs_pfvcolon_xlsx)

sessionInfo()
