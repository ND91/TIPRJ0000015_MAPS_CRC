#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop(paste0("Script needs 11 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
celltype_level <- args[2]
stackedbarplot_pbmc_pdf <- args[3]
stackedbarplot_pf_pdf <- args[4]
stackedbarplot_pbmc_pf_pdf <- args[5]
stackedbarplot_pbmc_myeloid_pdf <- args[6]
stackedbarplot_pf_myeloid_pdf <- args[7]
stackedbarplot_pbmc_pf_myeloid_pdf <- args[8]
stackedbarplot_pbmc_t_pdf <- args[9]
stackedbarplot_pf_t_pdf <- args[10]
stackedbarplot_pbmc_pf_t_pdf <- args[11]

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

seuratObject <- readRDS(seurat_rds_path)

# All

proportions_df <- seuratObject@meta.data %>%
  dplyr::group_by(UQ(rlang::sym(celltype_level)), SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(celltype_level))), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))

## PBMC

stackedbarplot_pbmc <- proportions_df %>%
  dplyr::filter(Tissue == "PBMC") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to CD45+ in PBMC",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## PF

stackedbarplot_pf <- proportions_df %>%
  dplyr::filter(Tissue == "PF") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to CD45+ in PF",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Together

stackedbarplot_pbmc_pf <- proportions_df %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to CD45+ in PBMCs",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~Tissue, nrow = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(width = 10, height = 5, file = stackedbarplot_pbmc_pdf)
print(stackedbarplot_pbmc)
dev.off()

pdf(width = 10, height = 5, file = stackedbarplot_pf_pdf)
print(stackedbarplot_pf)
dev.off()

pdf(width = 10, height = 5, file = stackedbarplot_pbmc_pf_pdf)
print(stackedbarplot_pbmc_pf)
dev.off()

# Myeloid

proportions_myeloid_df <- seuratObject@meta.data %>%
  dplyr::filter(manual_l1 == "Myeloid") %>%
  dplyr::group_by(UQ(rlang::sym(celltype_level)), SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(celltype_level))), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))

## PBMC

stackedbarplot_myeloid_pbmc <- proportions_myeloid_df %>%
  dplyr::filter(Tissue == "PBMC") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to myeloid in PBMC",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## PF

stackedbarplot_myeloid_pf <- proportions_myeloid_df %>%
  dplyr::filter(Tissue == "PF") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to myeloid in PF",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Together

stackedbarplot_myeloid_pbmc_pf <- proportions_myeloid_df %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to myeloid",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~Tissue, nrow = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(width = 7.5, height = 5, file = stackedbarplot_pbmc_myeloid_pdf)
print(stackedbarplot_myeloid_pbmc)
dev.off()

pdf(width = 7.5, height = 5, file = stackedbarplot_pf_myeloid_pdf)
print(stackedbarplot_myeloid_pf)
dev.off()

pdf(width = 7.5, height = 5, file = stackedbarplot_pbmc_pf_myeloid_pdf)
print(stackedbarplot_myeloid_pbmc_pf)
dev.off()

# T

proportions_t_df <- seuratObject@meta.data %>%
  dplyr::filter(manual_l1 == "T") %>%
  dplyr::group_by(UQ(rlang::sym(celltype_level)), SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(celltype_level))), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))

## PBMC

stackedbarplot_t_pbmc <- proportions_t_df %>%
  dplyr::filter(Tissue == "PBMC") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to T in PBMC",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## PF

stackedbarplot_t_pf <- proportions_t_df %>%
  dplyr::filter(Tissue == "PF") %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to T in PF",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Together

stackedbarplot_t_pbmc_pf <- proportions_t_df %>%
  ggplot(aes(x = Donor, y = Ncellprop)) +
  geom_bar(position="stack", stat="identity", aes(fill = UQ(rlang::sym(celltype_level)))) +
  labs(subtitle = "Proportion relative to T",
       y = "Proportion") +
  scale_fill_manual(values = col_vector) +
  facet_wrap(~Tissue, nrow = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(width = 10, height = 5, file = stackedbarplot_pbmc_t_pdf)
print(stackedbarplot_t_pbmc)
dev.off()

pdf(width = 10, height = 5, file = stackedbarplot_pf_t_pdf)
print(stackedbarplot_t_pf)
dev.off()

pdf(width = 10, height = 5, file = stackedbarplot_pbmc_pf_t_pdf)
print(stackedbarplot_t_pbmc_pf)
dev.off()


sessionInfo()