#!/usr/bin/env Rscript
# The goal of this script is to extract the CDR3s from the T lineage from the paired TX patients only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

airr_csv <- args[1]
crc_pbmc_pf_tx_t_paired_cdr3_airr_csv <- args[2]

airr <- read.csv(airr_csv)

txdonors <- airr %>%
  dplyr::filter(Tissue == "TX") %>%
  dplyr::pull(Donor) %>%
  unique()

airr_paired_tx_donor <- airr %>%
  dplyr::filter(Donor %in% txdonors)

# Save data
write.csv(airr_paired_tx_donor, crc_pbmc_pf_tx_t_paired_cdr3_airr_csv)

sessionInfo()