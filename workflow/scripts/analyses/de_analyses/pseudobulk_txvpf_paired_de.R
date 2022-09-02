#!/usr/bin/env Rscript
# This script will perform pseudobulk DE analysis using DESeq2.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
celltype_column_name <- args[2]
celltype_name <- args[3]
deg_csv_path <- args[4]
deseq2_object_path <- args[5]

seuratObject <- readRDS(seurat_rds_path)

if(!celltype_name %in% seuratObject@meta.data[,celltype_column_name]) stop(paste0(celltype_name, " cannot be found in SeuratObject metadata column: ", celltype_column_name))

sampleinfo <- unique(seuratObject@meta.data[,c("Donor", "Tissue")])
rownames(sampleinfo) <- paste0(sampleinfo$Donor, "_", sampleinfo$Tissue)

seuratObject <- seuratObject[,seuratObject@meta.data[,celltype_column_name] == celltype_name]
seuratObject@meta.data$DonorTissue <- paste0(seuratObject@meta.data$Donor, "_", seuratObject@meta.data$Tissue)

counts_mat <- GetAssayData(seuratObject, slot = "counts")

counts_sample <- counts_mat %*% model.matrix(~0+DonorTissue, data = seuratObject@meta.data)
colnames(counts_sample) <- gsub("DonorTissue", "", colnames(counts_sample))

counts_sample <- counts_sample[which(Matrix::rowSums(counts_sample) != 0),]

dds <- DESeqDataSetFromMatrix(countData = counts_sample,
                              colData = sampleinfo[colnames(counts_sample),],
                              design = ~Tissue)
dds <- DESeq(dds)

degs <- results(dds, contrast = c("Tissue", "TX", "PF"), independentFiltering = T)
degs <- degs[order(degs$pvalue),]
degs <- degs[!is.na(degs$padj),]
degs$gene <- rownames(degs)

write.csv(degs, deg_csv_path)
saveRDS(dds, deseq2_object_path, compress = "gzip")
