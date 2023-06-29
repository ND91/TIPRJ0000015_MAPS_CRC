#!/usr/bin/env Rscript

# The goal of this script is to create a stacked barplot of the cellular composition at l2 for for PBMC, PF, colon (public), and liver (public).

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

seurat_pbmc_pf_rds <- args[1]
seurat_liver_rds <- args[2]
seurat_colon_rds <- args[3]
marker_genes_xlsx <- args[4]
stackedbarplot_tissue_pdf <- args[5]
stackedbarplot_tissue_colonwoB_pdf <- args[6]
stackedbarplot_section_pdf <- args[7]

seuratObject_pbmc_pf <- readRDS(seurat_pbmc_pf_rds)
seuratObject_liver <- readRDS(seurat_liver_rds)
seuratObject_colon <- readRDS(seurat_colon_rds)
marker_genes <- readxl::read_excel(marker_genes_xlsx) %>%
  dplyr::filter(level == "manual_l2",
                !is.na(color)) %>%
  dplyr::select(celltype, color) %>%
  unique()

celltype_l2_color <- marker_genes$color
names(celltype_l2_color) <- marker_genes$celltype

head(seuratObject_liver@meta.data)

liver_metadata <- seuratObject_liver@meta.data %>%
  dplyr::rename(manual_l2 = cell_type_MAPS_L2,
                Donor = patient) %>%
  dplyr::mutate(SampleID = paste0(Donor, "_Liver")) %>%
  dplyr::mutate(Tissue = "Liver",
                Section = "Liver",
                Study = "[liver study]") %>%
  dplyr::select(Tissue, Section, manual_l2, Donor, SampleID, Study)
colon_metadata <- seuratObject_colon@meta.data %>%
  dplyr::rename(manual_l2 = cell_type_MAPS_L2,
                Donor = donor) %>%
  dplyr::mutate(Study = "[colon study]",
                Tissue = "Colon",
                Section = region,
                SampleID = paste0(Donor, "_", Tissue)) %>%
  dplyr::select(Tissue, Section, manual_l2, Donor, SampleID, Study)
pf_pbmc_metadata <- seuratObject_pbmc_pf@meta.data %>%
  dplyr::mutate(Study = "Current study",
                Section = Tissue) %>%
  dplyr::select(Tissue, Section, manual_l2, Donor, SampleID, Study)

combined_metadata <- rbind(liver_metadata, colon_metadata, pf_pbmc_metadata)

combined_metadata <- combined_metadata %>%
  dplyr::left_join(seuratObject_pbmc_pf@meta.data %>%
                     dplyr::select(manual_l1, manual_l2) %>%
                     unique())

# Per tissue

proportions_tissue_df <- combined_metadata %>%
  dplyr::group_by(manual_l2, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Tissue) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver")))

stackedbarplot_tissue_ggobj <- proportions_tissue_df %>%
  ggplot(aes(x = Tissue, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = manual_l2)) +
  labs(subtitle = "Proportion relative to CD45+",
       y = "Proportion") +
  scale_fill_manual(values = celltype_l2_color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 3, height = 5, file = stackedbarplot_tissue_pdf)
print(stackedbarplot_tissue_ggobj)
dev.off()

# Per tissue (Colon without B)

proportions_tissue_colonwoB_df <- combined_metadata %>%
  dplyr::filter(!(Tissue == "Colon" & manual_l2 == "Plasma B")) %>%
  dplyr::group_by(manual_l2, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Tissue) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver")))

stackedbarplot_tissue_colonwoB_ggobj <- proportions_tissue_colonwoB_df %>%
  ggplot(aes(x = Tissue, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = manual_l2)) +
  labs(subtitle = "Proportion relative to CD45+",
       y = "Proportion") +
  scale_fill_manual(values = celltype_l2_color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 3, height = 5, file = stackedbarplot_tissue_colonwoB_pdf)
print(stackedbarplot_tissue_colonwoB_ggobj)
dev.off()

# Per section

proportions_section_df <- combined_metadata %>%
  dplyr::group_by(manual_l2, Section, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Section), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Section) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Section = factor(Section, levels = c("PBMC", "PF", "Caecum", "Transverse", "Sigmoid", "Liver")))

stackedbarplot_section_ggobj <- proportions_section_df %>%
  ggplot(aes(x = Section, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = manual_l2)) +
  labs(subtitle = "Proportion relative to CD45+",
       y = "Proportion") +
  scale_fill_manual(values = celltype_l2_color) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 5, height = 5, file = stackedbarplot_section_pdf)
print(stackedbarplot_section_ggobj)
dev.off()

# Myeloid (THIS ONE STILL NEEDS TO AWAIT THE L3)

# proportions_myeloid_df <- combined_metadata %>%
#   dplyr::filter(manual_l1 == "Myeloid") %>%
#   dplyr::group_by(manual_l2, Tissue, .drop = T) %>% 
#   dplyr::summarise(Ncells = n()) %>%
#   dplyr::ungroup() %>%
#   tidyr::complete(tidyr::nesting(Tissue), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
#   dplyr::group_by(Tissue) %>%
#   dplyr::mutate(Nsample = sum(Ncells),
#                 Ncellprop = Ncells/Nsample,
#                 Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop))
# 
# stackedbarplot_myeloid_ggobj <- proportions_myeloid_df %>%
#   ggplot(aes(x = Tissue, y = Ncellprop)) +
#   geom_bar(position="stack", stat="identity", aes(fill = manual_l2)) +
#   labs(subtitle = "Proportion relative to myeloid",
#        y = "Proportion") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_blank(),
#         legend.pos = "bottom",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# pdf(width = 7.5, height = 5, file = stackedbarplot_pbmc_pf_myeloid_pdf)
# print(stackedbarplot_myeloid_ggobj)
# dev.off()

sessionInfo()