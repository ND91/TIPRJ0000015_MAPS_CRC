#!/usr/bin/env Rscript

# This script will perform marker gene expression analyses comparing the different macrophages at l3 between CRC PM+ with HC.

require(dplyr)
require(ggplot2)
require(ggrastr)
require(ggrepel)
require(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

deseq2_list_rds <- args[1]
volcanoplot_pdf <- args[2]

deseq2_list <- readRDS(deseq2_list_rds)

degs_df <- do.call(rbind, lapply(names(deseq2_list), function(celltype){
  data.frame(deseq2_list[[celltype]]$degs, celltype = celltype)
}))

degs_sig <- degs_df %>%
  dplyr::filter(padj<0.05)

genes_of_interest <- c("SLC2A3", "HBEGF", "HIF1A", "THBS1", "VEGFA", "IER2", "EREG", "IL10")

volcanoplot_ggobj <- degs_df %>%
  dplyr::mutate(Significance = ifelse(padj<0.05, "Significant", "NS"),
                celltype = factor(celltype, levels = c("M1-like", "M2-like", "M1/M2-like")),
                label = ifelse(gene %in% genes_of_interest, gene, NA)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_point_rast(aes(col = Significance)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(gene %in% genes_of_interest), col = "red") +
  geom_label_repel(aes(label = label)) +
  scale_color_manual(values = c("Significant" = "#000000", "NS" = "#d3d3d3")) +
  labs(title = "PM+CRC relative to HC",
       y = bquote('-'~log[10]~'(p-value)'),
       x = bquote('-'~log[2]~'(fold change)')) +
  facet_wrap(~celltype, nrow = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

pdf(width = 15, height = 5, file = volcanoplot_pdf)
print(volcanoplot_ggobj)
dev.off()

sessionInfo()