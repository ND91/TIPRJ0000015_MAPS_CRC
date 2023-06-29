#!/usr/bin/env Rscript

# This script will perform the differential expression analyses comparing CRC with HC for the desired tissue.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

seurat_rds <- args[1]
seuratDE_r <- args[2]
tissue <- args[3]
celltype_level <- args[4]
degs_list_rds <- args[5]
degs_xlsx <- args[6]

source(seuratDE_r)

seuratObject <- readRDS(seurat_rds)

group_celltype <- unique(seuratObject@meta.data[,c("Group", celltype_level)])
celltypes_present <- names(which(colSums(table(group_celltype$Group, group_celltype[,celltype_level]) != 0*1) == 2)) #Identify which celltypes are present in both tissues
seuratObject_filtered <- seuratObject[,seuratObject@meta.data[,celltype_level] %in% celltypes_present]
seuratObject_filtered <- seuratObject_filtered[,seuratObject_filtered@meta.data$Tissue == tissue]

pb_sample_metadata <- unique(seuratObject_filtered@meta.data[,c("SampleID", "Donor", "Group", "Sex", "Age", "Tissue")]) %>%
  dplyr::mutate(Group = ifelse(Group == "CRC+", "CRCp", Group),
                Group = ifelse(Group == "CRC-", "CRCn", Group))
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

degs_crcvhc <- seuratDE(seuratobj = seuratObject_filtered, 
                        cellsampleID = "SampleID", 
                        cellclusterID = celltype_level, 
                        sampleinfo = pb_sample_metadata, 
                        design = "~Group", 
                        contrast = c("Group", "CRCp", "HC"))

degs_crcvhc <- degs_crcvhc[lapply(degs_crcvhc, length) != 0]

saveRDS(degs_crcvhc, degs_list_rds, compress = "gzip")

write.xlsx(lapply(degs_crcvhc, function(celltype){
  celltype$degs
}), file = degs_xlsx)

sessionInfo()
