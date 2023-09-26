#!/usr/bin/env Rscript

# This script will perform the differential Nichenet analysis.

require(Seurat)
require(dplyr)
require(muscat)
require(openxlsx)
require(SingleCellExperiment)
require(nichenetr)
require(ggplot2)
require(viridis)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

seurat_rds <- args[1]
lr_network_rds <- args[2]
degs_crcvhc_PF_manual_l3_rds <- args[3]
celltype_markers_xlsx <- args[4]
ligand_target_matrix_rds <- args[5]
nichenet_output_rds <- args[6]

seuratObject <- readRDS(seurat_rds)
Idents(seuratObject) <- "manual_l3"
degs_crcvhc_PF_manual_l3 <- readRDS(degs_crcvhc_PF_manual_l3_rds)
lr_network <- readRDS(lr_network_rds)
ligand_target_matrix <- readRDS(ligand_target_matrix_rds)

degs_crcvhc_PF_manual_l3_sig <- lapply(degs_crcvhc_PF_manual_l3, function(gene){
  gene$degs %>%
    data.frame() %>%
    dplyr::filter(padj<0.05)
})

degs_crcvhc_PF_macrophages_vcan <- degs_crcvhc_PF_manual_l3_sig$`Macrophages SPP1+` %>%
  data.frame() %>%
  dplyr::filter(gene %in% rownames(ligand_target_matrix))

expressed_genes_macrophages_vcan <- get_expressed_genes("Macrophages SPP1+", seuratObject, pct = 0.10, assay = "RNA")
background_expressed_genes_macrophages_vcan <- expressed_genes_macrophages_vcan %>% .[. %in% rownames(ligand_target_matrix)]

nonmacrophages_vcan_celltypes <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l3) %>%
  summarize(ncount = n()) %>%
  dplyr::arrange(ncount) %>%
  dplyr::filter(ncount > 10)

list_expressed_genes_nonmacrophages_vcan <- vector("list", length(nonmacrophages_vcan_celltypes))
names(list_expressed_genes_nonmacrophages_vcan) <- nonmacrophages_vcan_celltypes

for(i in 1:nrow(nonmacrophages_vcan_celltypes)){
  list_expressed_genes_nonmacrophages_vcan[[i]] <- get_expressed_genes(nonmacrophages_vcan_celltypes$manual_l3[i], seurat_obj = seuratObject, pct = 0.1, assay = "RNA")
}

expressed_genes_nonmacrophages_vcan <- list_expressed_genes_nonmacrophages_vcan %>% 
  unlist() %>% 
  unique()

ligands <- lr_network %>% 
  dplyr::pull(from) %>% 
  unique()
receptors <- lr_network %>% 
  dplyr::pull(to) %>% 
  unique()


# 1. Which cells result in the differential expression of CRC+-associated DEGs in Macrophages VCAN+

expressed_ligands <- intersect(ligands, expressed_genes_macrophages_vcan)
expressed_receptors <- intersect(receptors, expressed_genes_nonmacrophages_vcan)

potential_ligands <- lr_network %>% 
  dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  dplyr::pull(from) %>% 
  unique()

laa_nonmacrophages_vcan_to_macrophages_vcan <- predict_ligand_activities(geneset = degs_crcvhc_PF_macrophages_vcan$gene, 
                                                                         background_expressed_genes = background_expressed_genes_macrophages_vcan, 
                                                                         ligand_target_matrix = ligand_target_matrix, 
                                                                         potential_ligands = potential_ligands)

laa_nonmacrophages_vcan_to_macrophages_vcan_sub <- laa_nonmacrophages_vcan_to_macrophages_vcan %>% 
  dplyr::arrange(aupr) %>% 
  dplyr::mutate(test_ligand = factor(test_ligand, levels = test_ligand)) %>%
  dplyr::slice_tail(n = 30)

barplot_laa_nonmacrophages_vcan_to_macrophages_vcan_laa_ggplotobj <- ggplot(laa_nonmacrophages_vcan_to_macrophages_vcan_sub, aes(x = aupr, y = test_ligand)) +
  geom_bar(stat = "identity") +
  labs(y = "Ligand",
       x = "AUPR") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Dotplot

laa_nonmacrophages_vcan_to_macrophages_vcan_degs <- do.call(rbind, lapply(names(degs_crcvhc_PF_manual_l3), function(celltype){
  degs_crcvhc_PF_manual_l3[[celltype]]$degs %>%
    data.frame() %>%
    dplyr::filter(gene %in% laa_nonmacrophages_vcan_to_macrophages_vcan_sub$test_ligand) %>%
    dplyr::mutate(celltype = celltype)
}))

dotplot_laa_nonmacrophages_vcan_to_macrophages_vcan_degs_ggplotobj <- laa_nonmacrophages_vcan_to_macrophages_vcan_degs %>%
  dplyr::mutate(direction = ifelse(stat<0, "HC", "CRC+"),
                significance = ifelse(pvalue<0.05, "Significant", "NS"),
                gene = factor(gene, levels = laa_nonmacrophages_vcan_to_macrophages_vcan_sub$test_ligand),
                celltype = factor(celltype, levels = manual_l3_order)) %>%
  ggplot(aes(x = celltype, y = gene)) +
  geom_point(aes(size = -log10(pvalue), col = direction, alpha = significance)) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


laa_lr_weight <- lapply(laa_nonmacrophages_vcan_to_macrophages_vcan_sub$test_ligand, 
       get_weighted_ligand_target_links,
       geneset = degs_crcvhc_PF_macrophages_vcan$gene, 
       ligand_target_matrix = ligand_target_matrix, 
       n = 200) %>% 
  dplyr::bind_rows() %>% 
  tidyr::drop_na()

laa_lr_weight_ggplotobj <- ggplot(laa_lr_weight, aes(x = target, y = ligand, fill = weight)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_bw() +
  labs(x = "Receptor", y = "Ligand") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("hc_crcpmp_pf_nn_macrophages_vcan.pdf", width = 16, height = 7.5)
pdf("hc_crcpmp_pf_nn_macrophages_c1q.pdf", width = 22.5, height = 7.5)
pdf("hc_crcpmp_pf_nn_macrophages_vcan_c1q.pdf", width = 16, height = 7.5)
pdf("hc_crcpmp_pf_nn_macrophages_spp1.pdf", width = 10, height = 7.5)
#annotate_figure(ggarrange(barplot_laa_nonmacrophages_vcan_to_macrophages_vcan_laa_ggplotobj, laa_lr_weight_ggplotobj, dotplot_laa_nonmacrophages_vcan_to_macrophages_vcan_degs_ggplotobj, nrow = 1, ncol = 3, align = "hv", widths = c(0.25, 0.5, 0.5)), top = text_grob("Macrophages SPP1+"))
annotate_figure(ggarrange(barplot_laa_nonmacrophages_vcan_to_macrophages_vcan_laa_ggplotobj, dotplot_laa_nonmacrophages_vcan_to_macrophages_vcan_degs_ggplotobj, nrow = 1, ncol = 2, align = "hv", widths = c(0.5, 1)), top = text_grob("Macrophages SPP1+"))
dev.off()

saveRDS(multinichenet_macrophagevcanvsall, multinichenet_macrophagevcanvsall_rds, compress = "gzip")


sessionInfo()
