#!/usr/bin/env Rscript

# This script will perform marker gene expression analyses comparing the different macrophages at l3 between PM+CRC with HC.

require(Seurat)
require(dplyr)
require(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seuratObject_rds <- args[1]
functions_r <- args[2]
deseq2_list_rds <- args[3]
degs_xlsx <-args[4]

source(functions_r)

seuratObject <- readRDS(seuratObject_rds)

selected_cells <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "PF",
                manual_l2 == "Macrophages")
seuratObject_selected <- seuratObject[,seuratObject@meta.data$cellID %in% selected_cells$cellID]

#Identify which celltypes are present in both tissues
group_celltype <- unique(seuratObject_selected@meta.data[,c("Group", "manual_l3")])
celltypes_present <- names(which(colSums(table(group_celltype$Group, group_celltype[,"manual_l3"]) != 0*1) == 2)) 
seuratObject_filtered <- seuratObject_selected[,seuratObject_selected@meta.data[,"manual_l3"] %in% celltypes_present]

pb_sample_metadata <- unique(seuratObject_filtered@meta.data[,c("SampleID", "Donor", "Group", "Sex", "Age")]) %>%
  dplyr::mutate(Group_recoded = factor(ifelse(Group == "CRC+", "CRCp", "HC"), levels = c("HC", "CRCp")))
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

degs_crcpvhc <- seuratDE(seuratobj = seuratObject_filtered, 
                         cellsampleID = "SampleID", 
                         cellclusterID = "manual_l3", 
                         sampleinfo = pb_sample_metadata, 
                         design = "~Group_recoded", 
                         contrast = c("Group_recoded", "CRCp", "HC"))

degs_crcpvhc <- degs_crcpvhc[lapply(degs_crcpvhc, length) != 0]

saveRDS(degs_crcpvhc, deseq2_list_rds, compress = "gzip")

openxlsx::write.xlsx(lapply(degs_crcpvhc, function(celltype){
  celltype$degs
}), file = degs_xlsx)

sessionInfo()