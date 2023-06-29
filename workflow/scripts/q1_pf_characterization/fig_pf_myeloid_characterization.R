#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_myeloid_rds <- args[1]
totalseq_xlsx <- args[2]
celltype_markers_xlsx <- args[3]
celltype_level
markers_of_interest_xlsx <- args[4]

#scatterplot_cite_gex_pdf <- args[3]

seuratObject <- readRDS(seurat_myeloid_rds)
markers_of_interest <- readxl::read_excel(markers_of_interest_xlsx) %>%
  dplyr::filter(Lineage == "Myeloid")

overlapping_gene_protein <- markers_of_interest %>%
  dplyr::count(Name) %>%
  dplyr::filter(n > 1)
markers_of_interest_pex <- markers_of_interest %>% 
  dplyr::filter(Modality == "PEX",
                Name %in% overlapping_gene_protein$Name)
markers_of_interest_gex <- markers_of_interest %>% 
  dplyr::filter(Modality == "GEX",
                Name %in% overlapping_gene_protein$Name)

celltype_order <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seuratObject@meta.data$manual_l3)) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

seuratObject <- NormalizeData(seuratObject, assay = "CITE", normalization.method = "CLR")
seuratObject <- NormalizeData(seuratObject, assay = "RNA", normalization.method = "CLR")

# Heatmap

marker_gex <- GetAssayData(seuratObject, assay = "RNA")[markers_of_interest_gex$Feature,]

marker_gex_median <- data.frame(FeatureID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-FeatureID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::filter(nUMIs != 0) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = mean(log1p(nUMIs))) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Celltype = factor(Celltype, unique(celltype_order$celltype)),
                FeatureID = factor(FeatureID, unique(markers_of_interest_gex$Feature)),
                Modality = "GEX") %>%
  dplyr::left_join(markers_of_interest_gex %>%
                     dplyr::select(-Modality), by = c("FeatureID" = "Feature"))

marker_pex <- GetAssayData(seuratObject, assay = "CITE")[markers_of_interest_pex$Feature,]

marker_pex_median <- data.frame(FeatureID = rownames(marker_pex), marker_pex) %>%
  tidyr::pivot_longer(-FeatureID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::filter(nUMIs != 0) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = mean(log1p(nUMIs))) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Celltype = factor(Celltype, unique(celltype_order$celltype)),
                FeatureID = factor(FeatureID, unique(markers_of_interest_pex$Feature)),
                Modality = "PEX") %>%
  dplyr::left_join(markers_of_interest_pex %>%
                     dplyr::select(-Modality), by = c("FeatureID" = "Feature"))

marker_median <- rbind(marker_gex_median, marker_pex_median)

feature_metadata <- marker_median %>%
  dplyr::ungroup() %>%
  dplyr::select(FeatureID, Modality, Name) %>%
  unique()

marker_median_wide <- marker_median %>%
  dplyr::select(-c(Modality, Lineage, Name)) %>%
  tidyr::pivot_wider(names_from = Celltype, values_from = Median) %>%
  data.frame(row.names = .$FeatureID, .) %>%
  dplyr::select(-FeatureID)


ComplexHeatmap::Heatmap(marker_median_wide[markers_of_interest %>%
                                             dplyr::filter(Name %in% overlapping_gene_protein$Name) %>%
                                             dplyr::pull(Feature),], 
                        name = "Myeloid", 
                        row_split = feature_metadata$Modality, 
                        #row_labels = feature_metadata$Name,
                        color_space = PurpleAndYellow(1000), 
                        cluster_columns = F, 
                        cluster_rows = F)
