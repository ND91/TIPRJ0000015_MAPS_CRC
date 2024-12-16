#!/usr/bin/env Rscript

# This script will perform nichetnet analysis interrogating the interactions within the PM-CRC TX samples where who is communicating with the macrophages.

require(Seurat)
require(dplyr)
require(nichenetr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

crcpmp_tx_seuratobject_rds <- args[1]
crcpmp_tx_marker_l2_list_rds <- args[2]
crcpmp_tx_marker_l3_list_rds <- args[3]
lr_network_rds <- args[4]
ligand_target_matrix_rds <- args[5]
weighted_networks_rds <- args[6]

seuratObject <- readRDS(crcpmp_tx_seuratobject_rds)
crcpmp_tx_marker_l2_list <- readRDS(crcpmp_tx_marker_l2_list_rds)
crcpmp_tx_marker_l3_list <- readRDS(crcpmp_tx_marker_l3_list_rds)
lr_network <- readRDS(lr_network_rds)
ligand_target_matrix <- readRDS(ligand_target_matrix_rds)
weighted_networks <- readRDS(weighted_networks_rds)

# seuratObject@meta.data$Celltype <- ifelse(seuratObject@meta.data$manual_l2 == "CD8T", "CD8T", seuratObject@meta.data$manual_l3)
seuratObject@meta.data$Celltype <- seuratObject@meta.data$manual_l2

Idents(seuratObject) <- "Celltype"

# crcpmp_tx_marker_list$CD8T$gene

expressed_genes_receiver <- get_expressed_genes(celltype_oi = "Macrophages", seuratObject, pct = 0.05)

# Define ligands from the receptors
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% 
  filter(to %in% expressed_receptors) %>% 
  pull(from) %>% 
  unique()

# 1. Agnostic sender

agnostic_sender_celltypes <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l2, manual_l3) %>%
  dplyr::summarize(cells = n()) %>%
  dplyr::filter(cells > 10,
                #manual_l2 != "CD8T"
                # manual_l3 != "Macrophages SPP1+"
                ) %>%
  dplyr::pull(manual_l2)

list_expressed_genes_agnostic_sender <- lapply(agnostic_sender_celltypes, get_expressed_genes, seuratObject, 0.05)
expressed_genes_agnostic_sender <- list_expressed_genes_agnostic_sender %>% 
  unlist() %>% 
  unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_agnostic_sender) 

geneset_oi <- crcpmp_tx_marker_l2_list$`Macrophages` %>% 
  dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25,
                gene %in% rownames(ligand_target_matrix))

background_expressed_genes <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))

ligand_activities_agnostic <- predict_ligand_activities(geneset = geneset_oi$gene,
                                                        background_expressed_genes = background_expressed_genes,
                                                        ligand_target_matrix = ligand_target_matrix,
                                                        potential_ligands = potential_ligands_focused)

ligand_activities_agnostic <- ligand_activities_agnostic %>% 
  dplyr::arrange(-aupr_corrected) %>% 
  dplyr::mutate(rank = rank(desc(aupr_corrected)))

ggplot(ligand_activities_agnostic, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities_agnostic %>% 
                                  dplyr::top_n(20, aupr_corrected) %>% 
                                  dplyr::pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_bw()

best_upstream_ligands <- ligand_activities_agnostic %>% 
  dplyr::top_n(20, aupr_corrected) %>% 
  dplyr::arrange(-aupr_corrected) %>% 
  pull(test_ligand)

vis_ligand_aupr <- ligand_activities_agnostic %>% 
  dplyr::filter(test_ligand %in% best_upstream_ligands) %>%
  tibble::column_to_rownames("test_ligand") %>% 
  dplyr::select(aupr_corrected) %>% 
  dplyr::arrange(aupr_corrected) %>% 
  as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))


active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi$gene,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% 
  tidyr::drop_na()


active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

## Ligands

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

## Receptors

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


p_dotplot <- DotPlot(seuratObject,
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# 2. Sender focused

ligand_activities_focused <- ligand_activities_agnostic %>% 
  dplyr::filter(test_ligand %in% potential_ligands_focused)

best_upstream_ligands <- ligand_activities_focused %>% 
  dplyr::top_n(30, aupr_corrected) %>% 
  dplyr::arrange(-aupr_corrected) %>%
  dplyr::pull(test_ligand) %>% 
  unique()

ligand_aupr_matrix <- ligand_activities_focused %>% 
  dplyr::filter(test_ligand %in% best_upstream_ligands) %>%
  tibble::column_to_rownames("test_ligand") %>% 
  dplyr::select(aupr_corrected) %>% 
  dplyr::arrange(aupr_corrected)

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi$gene,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% 
  tidyr::drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.5) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

ggsave("macrophagesreceiver_ligand_receptor.pdf", width = 12.5, height = 8)

p_dotplot <- DotPlot(subset(seuratObject),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right") + 
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9,  angle = 90, hjust = 0))

ggsave("macrophagesreceiver_ligand_dotplot.pdf", width = 5, height = 7)

# 2. Pre-specified gene set

pos_reg_t_chemotaxis <- c("ADAM10","ADAM17","CCL21","CCL5","CCR2","CXCL13","OXSR1","S100A7","STK39","WNK1","WNT5A","XCL1","XCL2")

geneset_oi <- data.frame(gene = intersect(pos_reg_t_chemotaxis, rownames(ligand_target_matrix)))
