#!/usr/bin/env Rscript
# The goal of this script is to convert the h5seurat object into a RDS as well as create a reannotated version thereof that is to be used in the analyses.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

pbmc_reference_h5seurat <- args[1]
pbmc_reference_rds <- args[2]

refPBMC <- LoadH5Seurat(pbmc_reference_h5seurat)

# Reannotate for merging later on.

refPBMC@meta.data <- refPBMC@meta.data %>%
  dplyr::select(-c(lane, time)) %>%
  dplyr::mutate(
    celltype.l4 = case_when(
      celltype.l3 %in% c("ASDC_mDC") ~ "ASDC MDC",
      celltype.l3 %in% c("ASDC_pDC") ~ "ASDC PDC",
      celltype.l3 %in% c("CD14 Mono") ~ "Classical monocyte",
      celltype.l3 %in% c("CD16 Mono") ~ "Non-classical monocyte",
      celltype.l3 %in% c("CD4 Naive") ~ "CD4 naive",
      celltype.l3 %in% c("CD4 Proliferating") ~ "CD4 proliferating",
      celltype.l3 %in% c("CD4 TCM_1") ~ "CD4 TCM 1",
      celltype.l3 %in% c("CD4 TCM_2") ~ "CD4 TCM 2",
      celltype.l3 %in% c("CD4 TCM_3") ~ "CD4 TCM 3",
      celltype.l3 %in% c("CD4 TEM_1") ~ "CD4 TEM 1",
      celltype.l3 %in% c("CD4 TEM_2") ~ "CD4 TEM 2",
      celltype.l3 %in% c("CD4 TEM_3") ~ "CD4 TEM 3",
      celltype.l3 %in% c("CD4 TEM_4") ~ "CD4 TEM 4",
      celltype.l3 %in% c("CD8 Naive") ~ "CD8 naive 1",
      celltype.l3 %in% c("CD8 Naive_2") ~ "CD8 naive 2",
      celltype.l3 %in% c("CD8 Proliferating") ~ "CD8 proliferating",
      celltype.l3 %in% c("CD8 TCM_1") ~ "CD8 TCM 1",
      celltype.l3 %in% c("CD8 TCM_2") ~ "CD8 TCM 2",
      celltype.l3 %in% c("CD8 TCM_3") ~ "CD8 TCM 3",
      celltype.l3 %in% c("CD8 TEM_1") ~ "CD8 TEM 1",
      celltype.l3 %in% c("CD8 TEM_2") ~ "CD8 TEM 2",
      celltype.l3 %in% c("CD8 TEM_3") ~ "CD8 TEM 3",
      celltype.l3 %in% c("CD8 TEM_4") ~ "CD8 TEM 4",
      celltype.l3 %in% c("CD8 TEM_5") ~ "CD8 TEM 5",
      celltype.l3 %in% c("CD8 TEM_6") ~ "CD8 TEM 6",
      celltype.l3 %in% c("cDC1") ~ "CDC1",
      celltype.l3 %in% c("cDC2_1") ~ "CDC2 1",
      celltype.l3 %in% c("cDC2_2") ~ "CDC2 2",
      celltype.l3 %in% c("dnT_1") ~ "DNT 1",
      celltype.l3 %in% c("dnT_2") ~ "DNT 2",
      celltype.l3 %in% c("Doublet") ~ "Multiplet",
      celltype.l3 %in% c("Eryth") ~ "Erythroblast",
      celltype.l3 %in% c("gdT_1") ~ "GDT 1",
      celltype.l3 %in% c("gdT_2") ~ "GDT 2",
      celltype.l3 %in% c("gdT_3") ~ "GDT 3",
      celltype.l3 %in% c("gdT_4") ~ "GDT 4",
      celltype.l3 %in% c("NK Proliferating") ~ "NK proliferating",
      celltype.l3 %in% c("NK_1") ~ "NK 1",
      celltype.l3 %in% c("NK_2") ~ "NK 2",
      celltype.l3 %in% c("NK_3") ~ "NK 3",
      celltype.l3 %in% c("NK_4") ~ "NK 4",
      celltype.l3 %in% c("NK_CD56bright") ~ "NK CD56",
      celltype.l3 %in% c("pDC") ~ "PDC",
      celltype.l3 %in% c("Platelet") ~ "Platelet",
      celltype.l3 %in% c("Treg Memory") ~ "CD4 Treg memory",
      celltype.l3 %in% c("Treg Naive") ~ "CD4 Treg naive",
      TRUE ~ celltype.l3),
    celltype.l3 = case_when(
      celltype.l4 %in% c("ASDC MDC") ~ "MDC",
      celltype.l4 %in% c("ASDC PDC") ~ "PDC",
      celltype.l4 %in% c("CD4 TCM 1", "CD4 TCM 2", "CD4 TCM 3") ~ "CD4 TCM",
      celltype.l4 %in% c("CD4 TEM 1", "CD4 TEM 2", "CD4 TEM 3", "CD4 TEM 4") ~ "CD4 TEM",
      celltype.l4 %in% c("CD8 naive 1", "CD8 naive 2") ~ "CD8 naive",
      celltype.l4 %in% c("CD8 TCM 1", "CD8 TCM 2", "CD8 TCM 3") ~ "CD8 TCM",
      celltype.l4 %in% c("CD8 TEM 1", "CD8 TEM 2", "CD8 TEM 3", "CD8 TEM 4", "CD8 TEM 5", "CD8 TEM 6") ~ "CD8 TEM",
      celltype.l4 %in% c("DNT 1", "DNT 2") ~ "DNT",
      celltype.l4 %in% c("GDT 1", "GDT 2", "GDT 3", "GDT 4") ~ "GDT",
      celltype.l4 %in% c("NK 1", "NK 2", "NK 3", "NK 4", "NK CD56") ~ "NK", 
      celltype.l4 %in% c("CDC2 1", "CDC2 2") ~ "CDC2", 
      celltype.l4 %in% c("CD4 Treg memory", "CD4 Treg naive") ~ "CD4 Treg",
      TRUE ~ celltype.l4),
    celltype.l2 = case_when(
      celltype.l3 %in% c("B intermediate kappa", "B intermediate lambda") ~ "B intermediate",
      celltype.l3 %in% c("B memory kappa", "B memory lambda") ~ "B memory",
      celltype.l3 %in% c("B naive kappa", "B naive lambda") ~ "B naive",
      celltype.l3 %in% c("Plasma", "Plasmablast") ~ "Plasma B",
      celltype.l3 %in% c("Classical monocyte", "Non-classical monocyte") ~ "Monocyte",
      celltype.l3 %in% c("CD4 CTL", "CD4 naive", "CD4 proliferating", "CD4 TCM", "CD4 TEM", "CD4 Treg") ~ "CD4T",
      celltype.l3 %in% c("CD8 naive", "CD8 proliferating", "CD8 TCM", "CD8 TEM") ~ "CD8T",
      celltype.l3 %in% c("CDC1", "CDC2") ~ "CDC",
      celltype.l3 %in% c("DNT", "MAIT", "GDT") ~ "other T",
      celltype.l3 %in% c("NK proliferating", "NK") ~ "NK",
      TRUE ~ celltype.l3),
    celltype.l1 = case_when(
      celltype.l2 %in% c("CDC", "PDC", "Platelet") ~ "Myeloid",
      celltype.l2 %in% c("B intermediate", "B memory", "B naive", "Plasma B") ~ "B",
      celltype.l2 %in% c("Monocyte", "CDC", "MDC") ~ "Myeloid",
      celltype.l2 %in% c("CD4T", "CD8T", "other T") ~ "T",
      celltype.l2 %in% c("Multiplet") ~ "Multiplet",
      celltype.l2 %in% c("NK") ~ "NK",
      celltype.l2 %in% c("ILC") ~ "ILC",
      TRUE ~ celltype.l2),
    celltype.l0 = "Immune",
    manual.l1 = celltype.l1,
    manual.l2 = celltype.l2,
    manual.l3 = celltype.l3,
    manual.l4 = celltype.l4) %>%
  dplyr::rename(SampleID = donor) %>%
  dplyr::mutate(Organism = "Homo sapiens",
                Group = "Normal",
                Tissue = "PBMC",
                Study = "Hao et al. 2021 (DOI: 10.1016/j.cell.2021.04.048)")

rownames(refPBMC@meta.data) <- colnames(refPBMC)

saveRDS(refPBMC, pbmc_reference_rds, compress = "gzip")

sessionInfo()
