#!/usr/bin/env Rscript

# The goal of this script is to create a UpsetR of the overlap between clonotypes.

require(ggplot2)
require(dplyr)
require(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

airr_csv <- args[1]
patient <- args[2]
upset_pdf <- args[3]

airr <- read.csv(airr_csv)

airr_tcr_libsize <- airr %>%
    dplyr::group_by(SampleID, .drop = T) %>% 
    dplyr::summarise(libsize = sum(consensus_count))

airr <- airr %>%
  dplyr::left_join(airr_tcr_libsize, by = "SampleID") %>%
  dplyr::mutate(clonotype_prop = consensus_count/libsize)

if(patient != "all"){
  airr_subset <- airr %>%
    dplyr::filter(Donor == patient)
} else{
  airr_subset <- airr
}

# Upset

clonotype_binary <- table(airr_subset$junction_aa, airr_subset$Tissue)==0
clonotype_cmobj <- make_comb_mat(clonotype_binary)
upset_plotobj <- UpSet(clonotype_cmobj)

pdf(width = 5, height = 3.5, file = upset_pdf)
print(upset_plotobj)
dev.off()

# Circos

sessionInfo()