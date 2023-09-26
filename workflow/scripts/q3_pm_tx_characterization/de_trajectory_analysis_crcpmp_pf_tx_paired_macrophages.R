#!/usr/bin/env Rscript

# This script will perform the DE analysis of the TSCAN trajectories found on PF- and TX-determined macrophages.

require(Seurat)
require(dplyr)
require(scater)
require(TSCAN)
require(tradeSeq)
require(slingshot)
require(scran)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

sce_rds <- args[1]
tscan_rootcl1_rds <- args[2]
degs_pathstat_xlsx <- args[3]
degs_tradeseq_rds <- args[4]

sce <- readRDS(sce_rds)
tscan_rootcl1_rds <- readRDS(tscan_rootcl1_rds)

# pathStat
sce$tscan_pf <- pathStat(tscan_rootcl1)[,"8"] # PF-lineage
sce$tscan_tx <- pathStat(tscan_rootcl1)[,"3"] # TX-lineage

sce_macrophagesvcan_pf <- sce[,colData(sce)$seurat_clusters == "1" & colData(sce)$Tissue == "PF"]

degs_pseudo_macrophagesvcan_txlin <- testPseudotime(sce_macrophagesvcan_pf, df=1, pseudotime=sce_macrophagesvcan_pf$tscan_tx)
degs_pseudo_macrophagesvcan_txlin <- degs_pseudo_macrophagesvcan_txlin[order(degs_pseudo_macrophagesvcan_txlin$p.value),]

degs_pseudo_macrophagesvcan_pflin <- testPseudotime(sce_macrophagesvcan_pf, df=1, pseudotime=sce_macrophagesvcan_pf$tscan_pf)
degs_pseudo_macrophagesvcan_pflin <- degs_pseudo_macrophagesvcan_pflin[order(degs_pseudo_macrophagesvcan_pflin$p.value),]

degs_pathstat <- list(TX = degs_pseudo_macrophagesvcan_txlin,
                      PF = degs_pseudo_macrophagesvcan_pflin)

openxlsx::write.xlsx(lapply(degs_pathstat, function(lineage){
  lineage
}), file = degs_pathstat_xlsx)

# tradeSeq
tscan_rootcl1_nona <- pathStat(tscan_rootcl1)
not_on_path <- is.na(tscan_rootcl1_nona)
tscan_rootcl1_nona[not_on_path] <- 0
cell_weights <- !not_on_path
storage.mode(cell_weights) <- "numeric"

tscan_rootcl1_gam <- fitGAM(counts(sce), 
                            pseudotime=tscan_rootcl1_nona,
                            cellWeights=cell_weights)

degs_tradeseq <- patternTest(tscan_rootcl1_gam)
degs_tradeseq <- degs_tradeseq[order(degs_tradeseq$pvalue),]

saveRDS(degs_tradeseq, degs_tradeseq_rds)