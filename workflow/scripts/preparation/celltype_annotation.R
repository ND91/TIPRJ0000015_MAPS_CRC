#!/usr/bin/env Rscript
# The goal of this script is to annotate cells using the PBMC reference set provided by Hao et al. 2021 alongside some other annotations

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
reference_rds_path <- args[2]
seurat_annotated_rds_path <- args[3]
seurat_annotated_csv_path <- args[4]

seuratObject <- readRDS(seurat_rds_path)

## Add % mitochondrial reads
seuratObject[["percent_MT"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")

## Add cell-cycle scoring (only if some cell cycle genes are present)
cc_genes_present <- list(
  s_genes = cc.genes.updated.2019$s.genes[which(cc.genes.updated.2019$s.genes %in% rownames(seuratObject))],
  g2m_genes = cc.genes.updated.2019$g2m.genes[which(cc.genes.updated.2019$g2m.genes %in% rownames(seuratObject))]
)

# For some reason when performing CellCycleScoring on SCTransformed data, it fails in some samples. The error is the following:
# Error in `cut_number()`:
#   ! Insufficient data values to produce 24 bins.
# Run `rlang::last_error()` to see where the error occurred.
# Having had a look on the Github page of Seurat, this seems to occur when no 24 bins can be found. An odd solution is to normalize with NormalizeData instead. 
# As the eventual goal is to merge the datasets and then renormalize later on again with SCTransform, I will use NormalizeData as a side-step.

if (any(unlist(lapply(cc_genes_present, length)) != 0)) {
  seuratObject <- CellCycleScoring(
    object = seuratObject,
    s.features = cc_genes_present$s_genes,
    g2m.features = cc_genes_present$g2m_genes
  )
}

## Add references from Human CITE-seq PF Atlas dataset. 
reference_seuratObject <- readRDS(reference_rds_path)

ref_anchors <- FindTransferAnchors(
  reference = reference_seuratObject,
  query = seuratObject,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

predictions <- TransferData(anchorset = ref_anchors , 
                            reference = reference_seuratObject,
                            refdata = list(
                              celltype_l1 = "celltype_l1",
                              celltype_l2 = "celltype_l2",
                              celltype_l3 = "celltype_l3",
                              celltype_l4 = "celltype_l4")
                            )
  
seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(celltype_l1 = predictions$celltype_l1$predicted.id,
                celltype_l1_score = predictions$celltype_l1$prediction.score.max, 
                celltype_l2 = predictions$celltype_l2$predicted.id, 
                celltype_l2_score = predictions$celltype_l2$prediction.score.max, 
                celltype_l3 = predictions$celltype_l3$predicted.id, 
                celltype_l3_score = predictions$celltype_l3$prediction.score.max,
                celltype_l4 = predictions$celltype_l4$predicted.id, 
                celltype_l4_score = predictions$celltype_l4$prediction.score.max)
rownames(seuratObject@meta.data) <- colnames(seuratObject)

cell_metadata <- seuratObject@meta.data

# Renormalize with SCTransform as this is the suggested normalization technique. Facilitates integration later on as well.

seuratObject <- DietSeurat(seuratObject, 
                           counts = TRUE, 
                           data = TRUE, 
                           scale.data = FALSE)

seuratObject <- SCTransform(seuratObject, verbose = FALSE, conserve.memory = TRUE)

# Save data
saveRDS(seuratObject, seurat_annotated_rds_path, compress = "gzip")
data.table::fwrite(cell_metadata, seurat_annotated_csv_path)

sessionInfo()