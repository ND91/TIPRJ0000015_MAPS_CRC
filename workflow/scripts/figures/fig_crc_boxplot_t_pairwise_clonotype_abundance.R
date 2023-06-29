#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the overlaps in unique clonotypes between the different tissues.

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrastr))
suppressPackageStartupMessages(require(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

airr_csv <- args[1]
boxplot_pdf <- args[2]
boxplot_patanno_pdf <- args[3]

airr <- read.csv(airr_csv)

airr_tcr_libsize <- airr %>%
    dplyr::group_by(SampleID, .drop = T) %>% 
    dplyr::summarise(libsize = sum(consensus_count))

airr <- airr %>%
  dplyr::left_join(airr_tcr_libsize, by = "SampleID") %>%
  dplyr::mutate(clonotype_prop = consensus_count/libsize)

clonotype_binary <- table(airr$junction_aa, airr$Tissue)==0

ol_frequency <- do.call(rbind, lapply(unique(airr$Donor), function(patient){
  airr_subset <- airr %>%
    dplyr::filter(Donor %in% patient)
  
  clonotype_binary <- table(airr_subset$junction_aa, airr_subset$Tissue)==0
  data.frame(Donor = patient,
             TXnPF = length(which(rowSums(clonotype_binary[,c("TX", "PF")]) == 2)),
             TXnPBMC = length(which(rowSums(clonotype_binary[,c("TX", "PBMC")]) == 2)),
             PBMCnPF = length(which(rowSums(clonotype_binary[,c("PBMC", "PF")]) == 2)),
             PBMCnPFnTX = length(which(rowSums(clonotype_binary) == 3)))
})) 

boxplot_obj <- ol_frequency %>%
  tidyr::pivot_longer(-Donor, names_to = "Comparison", values_to = "Frequencies") %>%
  dplyr::mutate(Comparison = factor(Comparison, levels = c("PBMCnPFnTX", "PBMCnPF", "TXnPBMC", "TXnPF")),
                Frequencies = Frequencies) %>%
  ggplot(aes(x = Comparison, y = Frequencies)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() + 
  labs(title = "Overlapping T clonotypes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 5, height = 5, file = boxplot_pdf)
print(boxplot_obj)
dev.off()

boxplot_patanno_obj <- ol_frequency %>%
  tidyr::pivot_longer(-Donor, names_to = "Comparison", values_to = "Frequencies") %>%
  dplyr::mutate(Comparison = factor(Comparison, levels = c("PBMCnPFnTX", "PBMCnPF", "TXnPBMC", "TXnPF")),
                Frequencies = Frequencies) %>%
  ggplot(aes(x = Comparison, y = Frequencies)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(col = Donor)) + 
  geom_line(aes(group = Donor, col = Donor), alpha = 0.25) +
  labs(title = "Overlapping T clonotypes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 5, height = 5, file = boxplot_patanno_pdf)
print(boxplot_patanno_obj)
dev.off()

sessionInfo()