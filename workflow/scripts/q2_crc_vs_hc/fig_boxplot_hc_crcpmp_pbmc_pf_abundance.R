#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
comparison <- args[2]
boxplot_pbmc_pdf <- args[3]
boxplot_pf_pdf <- args[4]
boxplot_pbmc_pf_pdf <- args[5]

seuratObject <- readRDS(seurat_rds_path)

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

if(reference_level != "manual_l0"){
  proportions_df <- seuratObject@meta.data %>%
    dplyr::filter(Tissue != "TX") %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level)), SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level))), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID, UQ(rlang::sym(reference_level))) %>%
    dplyr::mutate(Nlrefsample = sum(Ncells),
                  Ncellprop = Ncells/Nlrefsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                  Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
    dplyr::left_join(unique(seuratObject@meta.data[,c("SampleID", "Group")]), by = "SampleID")
  
  # PBMC
  
  boxplot_pbmc <- proportions_df %>%
    dplyr::filter(Tissue == "PBMC") %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage in PBMC",
         y = "Proportion") +
    scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # PF
  
  boxplot_pf <- proportions_df %>%
    dplyr::filter(Tissue == "PF") %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage in PF",
         y = "Proportion") +
    scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Together
  
  boxplot_pbmc_pf <- proportions_df %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
    facet_grid(Tissue~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage",
         y = "Proportion") +
    scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
  print(boxplot_pbmc)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pf_pdf)
  print(boxplot_pf)
  dev.off()
  
  pdf(width = 10, height = 10, file = boxplot_pbmc_pf_pdf)
  print(boxplot_pbmc_pf)
  dev.off()
  
} else{
  proportions_df <- seuratObject@meta.data %>%
    dplyr::filter(Tissue != "TX") %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), SampleID, Tissue, Group, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, Group, SampleID), tidyr::nesting(UQ(rlang::sym(case_level))), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Nsample = sum(Ncells),
                  Ncellprop = Ncells/Nsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                  Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))
  
  if(case_level %in% c("manual_l2", "manual_l3", "manual_l4")){
    proportions_df <- proportions_df %>%
      dplyr::left_join(seuratObject@meta.data %>%
                         dplyr::select(manual_l1, UQ(rlang::sym(case_level))) %>%
                         unique())
    
    # PBMC
    
    boxplot_pbmc_facetted <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # PF
    
    boxplot_pf_facetted <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Together
    
    boxplot_pbmc_pf_facetted <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      facet_grid(Tissue~manual_l1, scales = "free", space='free') +
      scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
    print(boxplot_pbmc_facetted)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_pdf)
    print(boxplot_pf_facetted)
    dev.off()
    
    pdf(width = 10, height = 10, file = boxplot_pbmc_pf_pdf)
    print(boxplot_pbmc_pf_facetted)
    dev.off()
    
  } else{
    # PBMC
    
    boxplot_pbmc <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # PF
    
    boxplot_pf <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Together
    
    boxplot_pbmc_pf <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      facet_grid(Tissue~., scales = "free", space='free') +
      #facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
    print(boxplot_pbmc)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_pdf)
    print(boxplot_pf)
    dev.off()
    
    pdf(width = 10, height = 10, file = boxplot_pbmc_pf_pdf)
    print(boxplot_pbmc_pf)
    dev.off()
  }
}

sessionInfo()