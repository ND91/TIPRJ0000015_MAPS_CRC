args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

counts_raw_path <- args[1]
counts_filter_path <- args[2]
max_ambiant_path <- args[3]
ambiant_stats_path <- args[4]

# library(DropletUtils)
# library(scater)
# library(Matrix)
# library(data.table)

print(paste0(Sys.time(), "\treading raw counts"))
counts_raw <- DropletUtils::read10xCounts(counts_raw_path)
print(paste0(Sys.time(), "\treading filtered counts"))
counts_filter <- DropletUtils::read10xCounts(counts_filter_path)
print(paste0(Sys.time(), "\testimateAmbience"))
ambient_genes <- DropletUtils::estimateAmbience(SingleCellExperiment::counts(counts_raw),
  good.turing = FALSE, round = FALSE
)
rm(counts_raw)

print(paste0(Sys.time(), "\tmaximumAmbience"))
sce_ambient <- DropletUtils::maximumAmbience(SingleCellExperiment::counts(counts_filter),
  ambient_genes,
  mode = "proportion"
)

features10x <- data.table::fread(file.path(counts_raw_path, "features.tsv.gz"), sep = "\t", header = F)

print(paste0(Sys.time(), "\twriting table "))
data.table::fwrite(
  x = data.frame(
    ensemble_gene = rownames(sce_ambient),
    hgnc_symbol = features10x$V2[match(rownames(sce_ambient), features10x$V1)],
    contamination = Matrix::rowMeans(sce_ambient, na.rm = TRUE)
  ),
  file = ambiant_stats_path,
  row.names = F,
  compress = "gzip"
)

print(paste0(Sys.time(), "\tSaving sce_ambiant RDS"))
saveRDS(sce_ambient, max_ambiant_path, compress = "gzip")
print(paste0(Sys.time(), "\t printing SessionInfo"))
sessionInfo()
