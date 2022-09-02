args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

seurat_rds_path <- args[1]
plot_df_csv_path <- args[2]
fig_svg_path <- args[3]

seuratObject <- readRDS(seurat_rds_path)

#Extract manual_l2: Macrophages and CD4T
seuratObject <- seuratObject[,seuratObject@meta.data$manual_l2 %in% c("Macrophages", "CD4T")]

l3rl2_abundances <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l2, manual_l3, Tissue, Donor) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::group_by(manual_l2, Tissue, Donor) %>%
  dplyr::mutate(Nl1sample = sum(Ncells),
                Ncellprop = Ncells/Nl1sample) %>%
  dplyr::filter(manual_l3 %in% c("M2-like", "CD4 Treg memory"),
                Tissue %in% c("PF", "TX"))

write.csv(l3rl2_abundances, plot_df_csv_path)

l3rl2_abundances_plot <- ggplot(l3rl2_abundances, aes(x = Tissue, y = Ncellprop)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  facet_grid(.~manual_l3, scales = "free", space='free') +
  labs(title = "TX vs PF",
       subtitle = "Proportion relative to the parent population (l3 relative l2)",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

svglite(width = 15, height = 7.5, file = fig_svg_path, bg = "white")
print(l3rl2_abundances_plot)
dev.off()

sessionInfo()