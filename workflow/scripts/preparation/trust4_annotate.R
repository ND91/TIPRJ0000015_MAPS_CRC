#!/usr/bin/env Rscript
# The goal of this script is to reannotate the TRUST4 output with the run and cell ID.

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

trust4_airr <- args[1]
runid <- args[2]
trust4_reannotated_airr <- args[3]

cdr3 <- read.csv(trust4_airr, sep = "\t") %>%
  dplyr::mutate(RunID = runid,
                cellID = paste0("S", gsub("-", "_", RunID), "_", gsub("([A-Z]+-[0-9])_.+$", "\\1", sequence_id)))

write.csv(cdr3, trust4_reannotated_airr, row.names = F)

sessionInfo()