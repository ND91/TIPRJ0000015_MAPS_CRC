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

degs_list_rds <- args[1]
volcanoplot_de_pdf <- args[2]
boxplot_de_pdf <- args[3]

degs_list <- readRDS(degs_list_rds)

genes_of_interest <- c("C1QA", "SPP1", "MAF", "CD163", "VCAN")

# Volcanoplot

degs_list$degs %>%
  data.frame() %>%
  dplyr::mutate(label = ifelse(gene %in% genes_of_interest, gene, NA),
                Significance = ifelse(padj<0.05, "Significant", "NS")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_hline(yintercept = 0, col = "#BEBEBE") +
  geom_vline(xintercept = 0, col = "#BEBEBE") +
  geom_point_rast(aes(col = Significance)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(!is.na(label)),
                  col = "red") +
  geom_label_repel(aes(label = label)) +
  labs(title = "M1-like vs M2-like",
       y = bquote('-'~log[10]~'(p-value)'),
       x = bquote('-'~log[2]~'(fold change)')) +
  scale_color_manual(values = c("NS" = "#d3d3d3", "Significant" = "#000000")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

# Boxplots

rld <- rlog(degs_list$dds)

rld %>%
  assay() %>%
  data.frame(., GeneID = rownames(.)) %>%
  dplyr::filter(GeneID %in% genes_of_interest) %>%
  tidyr::pivot_longer(-GeneID, names_to = "Donor_sample", values_to = "exprs") %>%
  dplyr::left_join(data.frame(colData(rld), Donor_sample = colnames(rld)), by = "Donor_sample") %>%
  dplyr::left_join(data.frame(degs_list$degs), by = c("GeneID" = "gene")) %>%
  dplyr::mutate(label = paste0(GeneID, "\np-value = ", formatC(pvalue, format = "e", digits = 3))) %>%
  ggplot(aes(x = celltype, y = exprs)) +
  geom_boxplot(aes(fill = celltype, shape = celltype), outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(fill = celltype), size = 2, shape = 21, alpha = 0.5) +
  labs(y = log[2]~"(counts)") +
  facet_wrap(~label, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom")
  

sessionInfo()