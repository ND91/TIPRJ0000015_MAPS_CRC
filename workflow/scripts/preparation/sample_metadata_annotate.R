suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_RDS_path <- args[1]
sample_metadata_path <- args[2]
seurat_sample_metadata_annotated_RDS_path <- args[3]
seurat_sample_metadata_annotated_csv_path <- args[4]

seuratObject <- readRDS(seurat_RDS_path)

sample_metadata <- readxl::read_excel(sample_metadata_path)

## Add sample metadata and add extra columns introduced by HTOdemultiplex to facilitate merging later only
seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::left_join(sample_metadata, by = "RunID")

seuratObject@meta.data <- seuratObject@meta.data[,order(colnames(seuratObject@meta.data))]

rownames(seuratObject@meta.data) <- colnames(seuratObject)

DefaultAssay(seuratObject) <- "RNA"

saveRDS(seuratObject, seurat_sample_metadata_annotated_RDS_path, compress = "gzip")
data.table::fwrite(seuratObject@meta.data, seurat_sample_metadata_annotated_csv_path)

sessionInfo()
