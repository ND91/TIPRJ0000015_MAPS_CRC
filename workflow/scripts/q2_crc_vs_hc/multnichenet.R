#!/usr/bin/env Rscript

# This script will perform the differential Nichenet analysis.

devtools::install_github("saeyslab/nichenetr")

require(Seurat)
require(dplyr)
require(openxlsx)
require(nichenetr)
require(multinichenetr)
require(SingleCellExperiment)

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

lr_network <- readRDS(lr_network_rds)
ligand_target_matrix <- readRDS(ligand_target_matrix_rds)

seuratObject <- readRDS(seurat_rds)
sce <- as.SingleCellExperiment(seuratObject)
sce <- alias_to_symbol_SCE(sce, "human") %>% 
  makenames_SCE()

colData(sce)$Group_recoded <- ifelse(colData(sce)$Group == "CRC+", "CRCp", "HC")

abundance_expression_info = get_abundance_expression_info(sce = sce, 
                                                          sample_id = "Donor", 
                                                          group_id = "Group", 
                                                          celltype_id = "manual_l3", 
                                                          min_cells = 10, 
                                                          senders_oi = unique(colData(sce)$manual_l3), 
                                                          receivers_oi = "Macrophages VCAN+", 
                                                          lr_network = lr_network, 
                                                          batches = NA)

saveRDS(laa_noncoi_to_coi, nichenet_output_rds, compress = "gzip")

sessionInfo()
