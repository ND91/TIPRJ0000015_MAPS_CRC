#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the all CD45+ cells level for PBMC, PF, liver (public), and colon (public).

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrastr))
suppressPackageStartupMessages(require(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

seurat_pbmc_pf_rds <- args[1]
seurat_liver_rds <- args[2]
seurat_colon_rds <- args[3]
marker_genes_xlsx <- args[4]
boxplot_tissue_pdf <- args[5]
boxplot_section_pdf <- args[6]

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
  dplyr::group_by(manual_l2, Tissue, SampleID, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(SampleID, Tissue), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver")))

boxplot_tissue_ggobj <- proportions_tissue_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l2, -Ncellprop), y = Ncellprop, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to all CD45+ cells",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(width = 12.5, height = 5, file = boxplot_tissue_pdf)
print(boxplot_tissue_ggobj)
dev.off()

# Per section

proportions_section_df <- combined_metadata %>%
  dplyr::group_by(manual_l2, Section, SampleID, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(SampleID, Section), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Section = factor(Section, levels = c("PBMC", "PF", "Caecum", "Transverse", "Sigmoid", "Liver")))

boxplot_section_ggobj <- proportions_section_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l2, -Ncellprop), y = Ncellprop, col = Section)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to all CD45+ cells",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 12.5, height = 5, file = boxplot_section_pdf)
print(boxplot_section_ggobj)
dev.off()

sessionInfo()
