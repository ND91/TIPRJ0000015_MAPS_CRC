#!/usr/bin/env Rscript

# This script will perform the differential expression analyses comparing macrophages with monocytes.

require(Seurat)
require(dplyr)
require(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
deseq2_list_rds <- args[3]

seuratObject <- readRDS(seurat_rds)

seuratObject_monomac <- seuratObject[,seuratObject@meta.data$manual_l2 %in% c("Monocytes", "Macrophages")]

#### TO BE CONTINUED


pb_sample_metadata <- unique(seuratObject_macrophage@meta.data[,c("SampleID", "manual_l3", "Donor", "Sex", "Age")]) %>%
  dplyr::mutate(celltype = gsub("-like", "", manual_l3))
rownames(pb_sample_metadata) <- paste0(pb_sample_metadata$SampleID, "_", pb_sample_metadata$celltype)

counts_mat <- GetAssayData(seuratObject_macrophage, slot = "counts")

seuratObject_macrophage@meta.data$sourceID <- paste0(seuratObject_macrophage@meta.data$SampleID, "_", gsub("-like", "", seuratObject_macrophage@meta.data$manual_l3))

counts_sample <- counts_mat %*% model.matrix(~0+sourceID, data = seuratObject_macrophage@meta.data)
colnames(counts_sample) <- gsub("sourceID", "", colnames(counts_sample))

counts_sample <- counts_sample[which(Matrix::rowSums(counts_sample) != 0),]

dds <- DESeqDataSetFromMatrix(countData = counts_sample,
                              colData = pb_sample_metadata[colnames(counts_sample),],
                              design = ~celltype+Donor)
dds <- DESeq(dds)

degs <- results(dds, contrast = c("celltype", "M1", "M2"), independentFiltering = T)
degs <- degs[order(degs$pvalue),]
degs <- degs[!is.na(degs$padj),]
degs$gene <- rownames(degs)
degs_list <- list(degs = degs,
                  dds = dds)

saveRDS(degs_list, deseq2_list_rds, compress = "gzip")

sessionInfo()
