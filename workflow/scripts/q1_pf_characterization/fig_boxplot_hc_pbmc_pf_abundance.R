#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrastr))
suppressPackageStartupMessages(require(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0("Script needs 8 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
comparison <- args[2]
boxplot_pbmc_pdf <- args[3]
boxplot_pbmc_patanno_pdf <- args[4]
boxplot_pf_pdf <- args[5]
boxplot_pf_patanno_pdf <- args[6]
boxplot_pbmc_pf_pdf <- args[7]
boxplot_pbmc_pf_patanno_pdf <- args[8]

seuratObject <- readRDS(seurat_rds_path)

case_level <- paste0("manual_", strsplit(comparison, "r")[[1]][1])
reference_level <- paste0("manual_", strsplit(comparison, "r")[[1]][2])

if(reference_level != "manual_l0"){
  proportions_df <- seuratObject@meta.data %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level)), SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(case_level)), UQ(rlang::sym(reference_level))), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID, UQ(rlang::sym(reference_level))) %>%
    dplyr::mutate(Nlrefsample = sum(Ncells),
                  Ncellprop = Ncells/Nlrefsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                  Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))
  
  # PBMC
  
  boxplot_pbmc <- proportions_df %>%
    dplyr::filter(Tissue == "PBMC") %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000") +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage in PBMCs",
         y = "Proportion") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  boxplot_pbmc_patanno <- proportions_df %>%
    dplyr::filter(Tissue == "PBMC") %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000", aes(shape = Donor)) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage in PBMCs",
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
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD") +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage in PF",
         y = "Proportion") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  boxplot_pf_patanno <- proportions_df %>%
    dplyr::filter(Tissue == "PF") %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD", aes(shape = Donor)) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space="free") +
    labs(subtitle = "Proportion relative to lineage in PF",
         y = "Proportion") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Together
  
  boxplot_pbmc_pf <- proportions_df %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space="free") +
    labs(subtitle = "Proportion relative to lineage",
         y = "Proportion") +
    theme_bw() +
    scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  boxplot_pbmc_pf_patanno <- proportions_df %>%
    ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
    facet_grid(.~ eval(rlang::sym(reference_level)), scales = "free", space='free') +
    labs(subtitle = "Proportion relative to lineage",
         y = "Proportion") +
    theme_bw() +
    scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
  print(boxplot_pbmc)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pbmc_patanno_pdf)
  print(boxplot_pbmc_patanno)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pf_pdf)
  print(boxplot_pf)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pf_patanno_pdf)
  print(boxplot_pf_patanno)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pbmc_pf_pdf)
  print(boxplot_pbmc_pf)
  dev.off()
  
  pdf(width = 10, height = 5, file = boxplot_pbmc_pf_patanno_pdf)
  print(boxplot_pbmc_pf_patanno)
  dev.off()
  
} else{
  proportions_df <- seuratObject@meta.data %>%
    dplyr::group_by(UQ(rlang::sym(case_level)), SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(UQ(rlang::sym(case_level))), fill = list(Ncells = 0)) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Nsample = sum(Ncells),
                  Ncellprop = Ncells/Nsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                  Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))
  
  if(case_level %in% c("manual_l2", "manual_l3")){
    proportions_df <- proportions_df %>%
      dplyr::left_join(seuratObject@meta.data %>%
                         dplyr::select(manual_l1, UQ(rlang::sym(case_level))) %>%
                         unique())
    
    # PBMC
    
    boxplot_pbmc_facetted <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000") +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_facetted_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000", aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000") +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000", aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # PF
    
    boxplot_pf_facetted <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD") +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pf_facetted_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD", aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pf <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD") +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pf_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD", aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Together
    
    boxplot_pbmc_pf_facetted <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_pf_facetted_patanno <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      facet_grid(.~manual_l1, scales = "free", space='free') +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_pf <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_pf_patanno <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facetted.pdf", boxplot_pbmc_pdf))
    print(boxplot_pbmc_facetted)
    dev.off()
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facettedpatanno.pdf", boxplot_pbmc_pdf))
    print(boxplot_pbmc_facetted_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
    print(boxplot_pbmc)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_patanno_pdf)
    print(boxplot_pbmc_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facetted.pdf", boxplot_pf_pdf))
    print(boxplot_pf_facetted)
    dev.off()
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facetted_patanno.pdf", boxplot_pf_pdf))
    print(boxplot_pf_facetted_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_pdf)
    print(boxplot_pf)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_patanno_pdf)
    print(boxplot_pf_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facetted.pdf", boxplot_pbmc_pf_pdf))
    print(boxplot_pbmc_pf_facetted)
    dev.off()
    
    pdf(width = 10, height = 5, file = gsub("\\.pdf", "_facetted_patanno.pdf", boxplot_pbmc_pf_pdf))
    print(boxplot_pbmc_pf_facetted_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pf_pdf)
    print(boxplot_pbmc_pf)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pf_patanno_pdf)
    print(boxplot_pbmc_pf_patanno)
    dev.off()
      
  } else{
    # PBMC
    
    boxplot_pbmc <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000") +
      labs(subtitle = "Proportion relative to all PBMC-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PBMC") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#F30000") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#F30000", aes(shape = Donor)) +
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
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD") +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pf_patanno <- proportions_df %>%
      dplyr::filter(Tissue == "PF") %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = "#224FBD")) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, col = "#224FBD") +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, col = "#224FBD", aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all PF-derived CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Together
    
    boxplot_pbmc_pf <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    boxplot_pbmc_pf_patanno <- proportions_df %>%
      ggplot(aes(x = forcats::fct_reorder(UQ(rlang::sym(case_level)), -Ncellprop), y = Ncellprop, col = Tissue)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
      labs(subtitle = "Proportion relative to all CD45+ cells",
           y = "Proportion") +
      theme_bw() +
      scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.pos = "bottom",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pdf)
    print(boxplot_pbmc)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_patanno_pdf)
    print(boxplot_pbmc_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_pdf)
    print(boxplot_pf)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pf_patanno_pdf)
    print(boxplot_pf_patanno)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pf_pdf)
    print(boxplot_pbmc_pf)
    dev.off()
    
    pdf(width = 10, height = 5, file = boxplot_pbmc_pf_patanno_pdf)
    print(boxplot_pbmc_pf_patanno)
    dev.off()
  }
}

sessionInfo()