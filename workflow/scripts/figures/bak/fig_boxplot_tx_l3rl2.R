args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

seurat_rds_path <- args[1]
plot_df_csv_path <- args[2]
fig_svg_path <- args[3]

seuratObject <- readRDS(seurat_rds_path)

#Extract manual_l1: B, T, Myeloid, and NK/ILC only
seuratObject <- seuratObject[,seuratObject@meta.data$manual_l1 %in% c("B", "T", "Myeloid", "NK/ILC")]

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l1 = ifelse(manual_l2 == "CD4T", "CD4T", manual_l1),
                manual_l1 = ifelse(manual_l2 == "CD8T", "CD8T", manual_l1),
                manual_l1 = ifelse(manual_l2 == "other T", "other T", manual_l1))

tx_l3rl1_abundances <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l1, manual_l3, Donor) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::group_by(Donor, manual_l1) %>%
  dplyr::mutate(Nl1sample = sum(Ncells),
                Ncellprop = Ncells/Nl1sample)

write.csv(tx_l3rl1_abundances, plot_df_csv_path)

tx_l3rl1_abundances_plot <- ggplot(tx_l3rl1_abundances, aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  facet_grid(.~manual_l1, scales = "free", space='free') +
  labs(title = "Tumor",
       subtitle = "Proportion relative to the lineage",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

svglite(width = 30, height = 7.5, file = fig_svg_path, bg = "white")
print(tx_l3rl1_abundances_plot)
dev.off()

sessionInfo()