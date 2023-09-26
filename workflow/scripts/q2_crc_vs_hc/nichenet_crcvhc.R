#!/usr/bin/env Rscript

# This script will perform the differential Nichenet analysis.

devtools::install_github("saeyslab/nichenetr")

require(Seurat)
require(dplyr)
require(openxlsx)
require(nichenetr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

seurat_rds <- args[1]
celltype <- args[2]
lr_network_rds <- args[3]
ligand_target_matrix_rds <- args[4]
degs_crcvhc_PF_manual_l3_rds <- args[5]
nichenet_output_rds <- args[6]

seuratObject <- readRDS(seurat_rds)
Idents(seuratObject) <- "manual_l3"
degs_crcvhc_PF_manual_l3 <- readRDS(degs_crcvhc_PF_manual_l3_rds)
lr_network <- readRDS(lr_network_rds)
ligand_target_matrix <- readRDS(ligand_target_matrix_rds)

degs_crcvhc_PF_manual_l3_sig <- lapply(degs_crcvhc_PF_manual_l3, function(gene){
  gene$degs %>%
    data.frame() %>%
    dplyr::filter(padj<0.05)
})

sort(unlist(lapply(degs_crcvhc_PF_manual_l3_sig, nrow)))

degs_crcvhc_PF_coi <- degs_crcvhc_PF_manual_l3_sig[[celltype]] %>%
  data.frame() %>%
  dplyr::filter(gene %in% rownames(ligand_target_matrix))

expressed_genes_coi <- get_expressed_genes(celltype, seuratObject, pct = 0.10, assay = "RNA")
background_expressed_genes_coi <- expressed_genes_coi %>% .[. %in% rownames(ligand_target_matrix)]

noncoi_celltypes <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l3) %>%
  summarize(ncount = n()) %>%
  dplyr::arrange(ncount) %>%
  dplyr::filter(ncount > 10)

list_expressed_genes_noncoi <- vector("list", length(noncoi_celltypes))
names(list_expressed_genes_noncoi) <- noncoi_celltypes

for(i in 1:nrow(noncoi_celltypes)){
  list_expressed_genes_noncoi[[i]] <- get_expressed_genes(noncoi_celltypes$manual_l3[i], seurat_obj = seuratObject, pct = 0.1, assay = "RNA")
}

expressed_genes_noncoi <- list_expressed_genes_noncoi %>% 
  unlist() %>% 
  unique()

ligands <- lr_network %>% 
  dplyr::pull(from) %>% 
  unique()
receptors <- lr_network %>% 
  dplyr::pull(to) %>% 
  unique()

expressed_ligands <- intersect(ligands, expressed_genes_coi)
expressed_receptors <- intersect(receptors, expressed_genes_noncoi)

potential_ligands <- lr_network %>% 
  dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  dplyr::pull(from) %>% 
  unique()

laa_noncoi_to_coi <- predict_ligand_activities(geneset = degs_crcvhc_PF_coi$gene, 
                                               background_expressed_genes = background_expressed_genes_coi, 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               potential_ligands = potential_ligands)

saveRDS(laa_noncoi_to_coi, nichenet_output_rds, compress = "gzip")

sessionInfo()
