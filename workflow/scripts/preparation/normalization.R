#!/usr/bin/env Rscript
# The goal of this script is to import and normalize single cell experiment sample

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

count_path <- args[1]
seurat_rds_path <- args[2]
runid <- args[3]

cnts <- Read10X(count_path)
if (class(cnts) == "list") {
  expr_cnts <- cnts$`Gene Expression`
  
  nonzero_cells <- colnames(expr_cnts)[colSums(expr_cnts) != 0]
  expr_cnts <- expr_cnts[,colnames(expr_cnts) %in% nonzero_cells]
  
  abc_cnts <- cnts$`Antibody Capture`
  colnames(abc_cnts) <- gsub("-[0-9]$", "", colnames(abc_cnts))

  seuratObject <- CreateSeuratObject(counts = expr_cnts, min.cells = 3)
  seuratObject[["HTO"]] <- CreateAssayObject(counts = abc_cnts)

  seuratObject <- NormalizeData(seuratObject)
} else {
  nonzero_cells <- colnames(cnts)[colSums(cnts) != 0]
  cnts <- cnts[,colnames(cnts) %in% nonzero_cells]
  
  seuratObject <- NormalizeData(CreateSeuratObject(counts = cnts, min.cells = 3))
}

# Append the run ID into the cellbarcode to prevent any collisions later on.
runid_rfriendly <- paste0("S", gsub("-", "_", runid))
seuratObject@meta.data$RunID <- runid
seuratObject@meta.data$Study <- "Current study"
seuratObject <- RenameCells(object = seuratObject, add.cell.id = runid_rfriendly)
seuratObject@meta.data$CellID <- rownames(seuratObject@meta.data)

saveRDS(seuratObject, seurat_rds_path, compress = "gzip")

sessionInfo()
