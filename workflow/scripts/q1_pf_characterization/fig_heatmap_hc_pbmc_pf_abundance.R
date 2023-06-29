#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

library(Seurat)
library(ComplexHeatmap)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
comparison <- args[2]
heatmap_pbmc_pf_scaled_pdf <- args[3]
heatmap_pbmc_pf_unscaled_pdf <- args[4]

seuratObject <- readRDS(seurat_rds_path)

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

if(reference_level != "manual_l0"){
  proportions_long <- seuratObject@meta.data %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level)), SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level))), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID, UQ(rlang::sym(reference_level))) %>%
    dplyr::mutate(Nlrefsample = sum(Ncells),
                  Ncellprop = Ncells/Nlrefsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
    dplyr::ungroup()
  
  proportions_wide <- proportions_long %>%
    dplyr::select(-c(Tissue, Nlrefsample, Ncells, eval(reference_level))) %>%
    tidyr::pivot_wider(names_from = SampleID, values_from = Ncellprop) %>%
    tibble::column_to_rownames(var = case_level)
  
  proportions_wide_scaled <- t(apply(proportions_wide*100, 1, scale))
  colnames(proportions_wide_scaled) <- colnames(proportions_wide)
  
  heatmap_colannotation <- proportions_long %>%
    dplyr::select(SampleID, Tissue) %>%
    unique() %>%
    tibble::column_to_rownames(var = "SampleID")
  
  heatmap_rowannotation <- proportions_long %>%
    dplyr::select(UQ(rlang::sym(reference_level)), UQ(rlang::sym(case_level))) %>%
    unique() %>%
    dplyr::arrange(UQ(rlang::sym(reference_level)), UQ(rlang::sym(case_level))) 
  
  heatmap_plotobj_scaled <- Heatmap(proportions_wide_scaled[unlist(heatmap_rowannotation[,case_level]),], 
                                    name = "Proportion relative to lineage",
                                    show_row_dend = F,
                                    split = unlist(heatmap_rowannotation[,reference_level]),
                                    cluster_columns = F,
                                    col = colorRampPalette(c("blue", "white", "red"))(100),
                                    #rect_gp = gpar(col = "black", lwd = 1),
                                    top_annotation = HeatmapAnnotation(Tissue = heatmap_colannotation$Tissue,
                                                                       gp = gpar(col = "black"),
                                                                       col = list(Tissue = c("PF" = "#224FBD", "PBMC" = "#F30000")),
                                                                       annotation_name_side = "left"))
  
  heatmap_plotobj_unscaled <- Heatmap(proportions_wide[unlist(heatmap_rowannotation[,case_level]),], 
                                      name = "Proportion relative to lineage",
                                      show_row_dend = F,
                                      split = unlist(heatmap_rowannotation[,reference_level]),
                                      cluster_columns = F,
                                      col = colorRampPalette(c("blue", "white", "red"))(100),
                                      #rect_gp = gpar(col = "black", lwd = 1),
                                      top_annotation = HeatmapAnnotation(Tissue = heatmap_colannotation$Tissue,
                                                                         gp = gpar(col = "black"),
                                                                         col = list(Tissue = c("PF" = "#224FBD", "PBMC" = "#F30000")),
                                                                         annotation_name_side = "left"))
} else{
  proportions_long <- seuratObject@meta.data %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, SampleID), UQ(rlang::sym(case_level)), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Nlrefsample = sum(Ncells),
                  Ncellprop = Ncells/Nlrefsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
    dplyr::ungroup()
  
  proportions_wide <- proportions_long %>%
    dplyr::select(-c(Tissue, Nlrefsample, Ncells)) %>%
    tidyr::pivot_wider(names_from = SampleID, values_from = Ncellprop) %>%
    tibble::column_to_rownames(var = case_level)
  
  proportions_wide_scaled <- t(apply(proportions_wide*100, 1, scale))
  colnames(proportions_wide_scaled) <- colnames(proportions_wide)
  
  heatmap_colannotation <- proportions_long %>%
    dplyr::select(SampleID, Tissue) %>%
    unique() %>%
    tibble::column_to_rownames(var = "SampleID")
  
  heatmap_rowannotation <- proportions_long %>%
    dplyr::select(UQ(rlang::sym(case_level))) %>%
    unique() %>%
    dplyr::arrange(UQ(rlang::sym(case_level))) 
  
  heatmap_plotobj_scaled <- Heatmap(proportions_wide_scaled[unlist(heatmap_rowannotation[,case_level]),], 
                             name = "Proportion relative to CD45+",
                             show_row_dend = F,
                             cluster_columns = F,
                             col = colorRampPalette(c("blue", "white", "red"))(100),
                             #rect_gp = gpar(col = "black", lwd = 1),
                             top_annotation = HeatmapAnnotation(Tissue = heatmap_colannotation$Tissue,
                                                                gp = gpar(col = "black"),
                                                                col = list(Tissue = c("PF" = "#224FBD", "PBMC" = "#F30000")),
                                                                annotation_name_side = "left"))
  
  heatmap_plotobj_unscaled <- Heatmap(proportions_wide[unlist(heatmap_rowannotation[,case_level]),], 
                                      name = "Proportion relative to CD45+",
                                      show_row_dend = F,
                                      cluster_columns = F,
                                      col = colorRampPalette(c("blue", "white", "red"))(100),
                                      #rect_gp = gpar(col = "black", lwd = 1),
                                      top_annotation = HeatmapAnnotation(Tissue = heatmap_colannotation$Tissue,
                                                                         gp = gpar(col = "black"),
                                                                         col = list(Tissue = c("PF" = "#224FBD", "PBMC" = "#F30000")),
                                                                         annotation_name_side = "left"))
}

if(case_level == "manual_l1"){
  pdf(width = 7, height = 2.5, file = heatmap_pbmc_pf_scaled_pdf)
  print(heatmap_plotobj_scaled)
  dev.off()
  
  pdf(width = 7, height = 2.5, file = heatmap_pbmc_pf_unscaled_pdf)
  print(heatmap_plotobj_unscaled)
  dev.off()
} else if(case_level == "manual_l2"){
  pdf(width = 7, height = 5, file = heatmap_pbmc_pf_scaled_pdf)
  print(heatmap_plotobj_scaled)
  dev.off()
  
  pdf(width = 7, height = 5, file = heatmap_pbmc_pf_unscaled_pdf)
  print(heatmap_plotobj_unscaled)
  dev.off()
} else if(case_level == "manual_l3"){
  pdf(width = 7, height = 7, file = heatmap_pbmc_pf_scaled_pdf)
  print(heatmap_plotobj_scaled)
  dev.off()
  
  pdf(width = 7, height = 7, file = heatmap_pbmc_pf_unscaled_pdf)
  print(heatmap_plotobj_unscaled)
  dev.off()
}

sessionInfo()
