#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(broom)
require(tidyr)
require(purrr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
totalseq_xlsx <- args[2]
scatterplot_cite_gex_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)
totalseq <- readxl::read_excel(totalseq_xlsx) %>%
  data.frame() %>%
  dplyr::filter(!is.na(Ensembl)) %>%
  dplyr::mutate(hgnc = mapIds(x = org.Hs.eg.db,
                              keys = Ensembl,
                              column = "SYMBOL", 
                              keytype = "ENSEMBL", 
                              multiVals = "first"))

seuratObject <- NormalizeData(seuratObject, assay = "CITE", normalization.method = "CLR")

cite_gex <- data.frame(GetAssayData(seuratObject, assay = "RNA")[which(rownames(GetAssayData(seuratObject, assay = "RNA")) %in% totalseq$hgnc),]) %>%
  tibble::rownames_to_column(var = "Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::left_join(totalseq, by = c("Gene" = "hgnc"), relationship = "many-to-many") %>%
  dplyr::select(-c(Cell_description, Antibody_description, clone, Sequence, Include, DNA_ID, Epithelial, Endothelial, T, NK.ILC, B, Myeloid))
cite_pex <- GetAssayData(seuratObject, assay = "CITE") %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "TSID") %>%
  tidyr::pivot_longer(-TSID, names_to = "CB", values_to = "nExpr")  %>%
  dplyr::left_join(totalseq, by = c("TSID" = "Name"), relationship = "many-to-many")

cite_total <- cite_gex %>%
  dplyr::inner_join(cite_pex, by = c("CB", "Ensembl"), relationship = "many-to-many")

cite_cor <- cite_total %>%
  nest(data = -Antibody_description) %>%
  dplyr::mutate(cor_df = map(data, ~ cor.test(.x$nUMIs, .x$nExpr)),
                tidied = map(cor_df, tidy)
  ) %>%
  dplyr::select(-c(data, cor_df)) %>%
  tidyr::unnest(tidied) %>%
  dplyr::arrange(estimate)

pdf(file = scatterplot_cite_gex_pdf, width = 25, height = 37.5)
cite_total %>% 
  dplyr::left_join(cite_cor, by = "Antibody_description") %>%
  dplyr::mutate(label = paste0(Antibody_description, "\ncor = ", round(estimate, 2))) %>%
  ggplot(aes(x = nUMIs, y = nExpr)) +
  geom_point_rast() +
  labs(y = "Normalized protein expression",
       x = "Normalized gene expression") + 
  facet_wrap(~label, ncol = 10, scales = "free") +
  theme_bw()
dev.off()  

sessionInfo()