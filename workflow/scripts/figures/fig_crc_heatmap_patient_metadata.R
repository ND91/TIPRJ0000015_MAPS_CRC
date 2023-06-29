#!/usr/bin/env Rscript

# The goal of this script is to create a heatmap of the CRC patient characteristics.

require(ComplexHeatmap)
require(GetoptLong)
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

patient_metadata_xlsx <- args[1]
heatmap_patient_metadata_pdf <- args[2]
heatmap_patient_metadata_nocytof_pdf <- args[3]

patient_metadata <- readxl::read_excel(patient_metadata_xlsx) %>%
  dplyr::rename("PatientID" = "patient nr",
                "CastorID" = "castor ID",
                "Age" = "age at surgery",
                "Sex" = "sex",
                "Group" = "group",
                "MMRprim" = "MMR status primary tumor (MMRp/MMRd)",
                "MSprim" = "MS status primary tumor (MSS/MSI)",
                "MMRPM" = "MMR status PM (MMRp/MMRd)",
                "MSPM" = "MS status PM (MSS/MSI)",
                "Omics" = "entity",
                "Tumor_omics" = "Tx entity",
                "Cancer_stage" = "Stage",
                "T_stage" = "T",
                "N_stage" = "N",
                "M_stage" = "M",
                "PM" = "PM status",
                "MAPS_CRC" = "in study 'CRC MAPS'",
                "MAPS_GC" = "in study 'GC MAPS'",
                "PF_ATLAS" = "in study 'atlas'") %>%
  dplyr::mutate(PatientID = paste0("pt", PatientID),
                Group = gsub("(\\+|\\-)", "", Group),
                Group = gsub(" \\(HIPEC\\)", "", Group),
                Group = gsub(" ", "", Group),
                BMI = gsub("NA", NA, BMI),
                BMI = as.numeric(BMI),
                Sex = stringr::str_to_title(Sex),
                PM = stringr::str_to_title(PM),
                PM = gsub("Na", NA, PM),
                PM = gsub("Pm\\+", "Yes", PM),
                PM = gsub("Pm\\-", "No", PM),
                PM = gsub("Hc", NA, PM),
                MAPS_CRC = stringr::str_to_title(MAPS_CRC),
                MAPS_GC = stringr::str_to_title(MAPS_GC),
                PF_ATLAS = stringr::str_to_title(PF_ATLAS),
                HIPEC = stringr::str_to_title(HIPEC),
                HIPEC = gsub("Hc", NA, HIPEC),
                Omics = gsub("both", "CyTOF+scRNAseq", Omics),
                Omics = gsub("scRNA-seq", "scRNAseq", Omics),
                Omics = gsub("CYTOF", "CyTOF", Omics),
                MMRprim = gsub(" - request", "", MMRprim),
                MMRprim = gsub(" \\(methyl\\)", "", MMRprim),
                MMRprim = gsub("results in reffering hospital", "NA", MMRprim),
                MMRprim = gsub("HC", NA, MMRprim),
                MMRprim = gsub("NA", "Unknown", MMRprim),
                MMRprim = gsub("Unclear from sample", "Unknown", MMRprim),
                MMRprim = gsub("not performed - no request", "Unknown", MMRprim),
                MMRprim = gsub("no resection", "Unknown", MMRprim),
                MSprim = gsub("NA", "Unknown", MSprim),
                MSprim = gsub("Unclear from sample", "Unknown", MSprim),
                MSprim = gsub("not performed - no request", "Unknown", MSprim),
                MSprim = gsub("no resection", "Unknown", MSprim),
                MMRPM = gsub("NA", "Unknown", MMRPM),
                MMRPM = gsub("Unclear from sample", "Unknown", MMRPM),
                MMRPM = gsub("no PM", "No PM", MMRPM),
                MMRPM = gsub("not performed - no request", "Unknown", MMRPM),
                MSPM = gsub("NA", "Unknown", MSPM),
                MSPM = gsub("not performed - no request", "Unknown", MSPM),
                Cancer_stage = gsub("(HC|NA)", NA, Cancer_stage),
                T_stage = gsub("(HC|NA)", NA, T_stage),
                N_stage = gsub("(HC|NA)", NA, N_stage),
                M_stage = gsub("(HC|NA)", NA, M_stage),
                Tumor_omics = gsub("both", "*#", Tumor_omics),
                Tumor_omics = gsub("scRNAseq", "*", Tumor_omics),
                Tumor_omics = gsub("CyTOF", "#", Tumor_omics),
                Tumor_omics = gsub("no", "", Tumor_omics),
                label = ifelse(!is.na(Tumor_omics), paste0(PatientID, Tumor_omics), PatientID),
                Group_summarized = case_when(
                  Group == "CRC" & PM == "Yes" ~ "CRC+",
                  Group == "CRC" & PM == "No" ~ "CRC-",
                  Group == "GC" & PM == "Yes" ~ "GC+",
                  Group == "GC" & PM == "No" ~ "GC-",
                  T ~ Group
                )
  ) %>%
  dplyr::select(label, MAPS_CRC, MAPS_GC, PF_ATLAS, PatientID, Tumor_omics, Sex, Age, BMI, Group, Group_summarized, Omics, MMRprim, MSprim, MMRPM, MSPM, Cancer_stage, T_stage, N_stage, M_stage, PM, HIPEC) %>%
  dplyr::filter(Group_summarized == "CRC+" & MAPS_CRC == "Yes" | Group_summarized == "HC",
                Omics != "CyTOF")

rownames(patient_metadata) <- patient_metadata$PatientID  

patient_metadata <- patient_metadata[order(patient_metadata$PM, patient_metadata$HIPEC, patient_metadata$PatientID),]

patient_metadata_labels <- mapply(function(x, y) {
  ifelse(!is.na(y), qq("@{x}<sup>@{y}</sup>"), qq("@{x}"))
  #qq("@{x}<sup>@{y}</sup>")
}, x = patient_metadata$PatientID, y = patient_metadata$Tumor_omics)

patient_metadata_plot <- HeatmapAnnotation(patientID = anno_text(x = gt_render(patient_metadata_labels, 
                                                                               align_widths = TRUE), 
                                                                 just = "left", 
                                                                 location = 0.25),
                                           Sex = patient_metadata$Sex,
                                           Age = anno_barplot(patient_metadata$Age, bar_width = 1, gp = gpar(fill = "#444444")),
                                           BMI = anno_barplot(patient_metadata$BMI, bar_width = 1, gp = gpar(fill = "#a69eb0")),
                                           Group = patient_metadata$Group_summarized, 
                                           HIPEC = patient_metadata$HIPEC,
                                           `MS primary tumor` = patient_metadata$MSprim,
                                           `MS peritoneal metastasis` = patient_metadata$MSPM,
                                           # PM = patient_metadata$PM,
                                           `Pathological stage` = anno_simple(x = patient_metadata$Cancer_stage, 
                                                                              pch = patient_metadata$Cancer_stage, 
                                                                              col = c("I" = "#ffffff",
                                                                                      "II" = "#ffffff",
                                                                                      "III" = "#ffffff",
                                                                                      "IV" = "#ffffff"),
                                                                              na_col = "#ffffff",
                                                                              gp = gpar(col = "black")),
                                           Omics = patient_metadata$Omics,
                                           annotation_name_side = "left", 
                                           gp = gpar(col = "black"),
                                           border = T,
                                           na_col = "#ffffff",
                                           col = list(Sex = c("Male" = "#3b55cd", 
                                                              "Female" = "#cf3476"),
                                                      Group = c("HC" = "#91d1c2",
                                                                "CRC-" = "#7f2284",
                                                                "CRC+" = "#845422"),
                                                      Omics = c("CyTOF+scRNAseq" = "#009999",
                                                                "scRNAseq" = "#FFFAA0",
                                                                "CITEseq" = "#A7C7E7"),
                                                      `MS primary tumor` = c("HC" = "#ffffff",
                                                                             "MSI" = "#ff0000",
                                                                             "MSS" = "#0000ff",
                                                                             "Unknown" = "#a69eb0"),
                                                      `MS peritoneal metastasis` = c("MSI" = "#ff0000",
                                                                                     "MSS" = "#0000ff",
                                                                                     "Unknown" = "#a69eb0",
                                                                                     "No PM" = "#ffffff"),
                                                      HIPEC = c("Yes" = "#00BFC4",
                                                                "No" = "#7CAE00")
                                                      # PM = c("No" = "#054f42",
                                                      #        "Yes" = "#AB2020")
                                                      )
                                           )

pdf(width = 11, height = 7, file = heatmap_patient_metadata_pdf)
plot(patient_metadata_plot)
dev.off()

# No CyTOF-only

patient_metadata_nocytof <- patient_metadata %>%
  dplyr::filter(Omics %in% c("scRNAseq", "CyTOF+scRNAseq"))

patient_metadata_nocytof_labels <- mapply(function(x, y) {
  ifelse(!is.na(y), qq("@{x}<sup>@{y}</sup>"), qq("@{x}"))
  #qq("@{x}<sup>@{y}</sup>")
}, x = patient_metadata_nocytof$PatientID, y = patient_metadata_nocytof$Tumor_omics)

patient_metadata_nocytof_plot <- HeatmapAnnotation(patientID = anno_text(x = gt_render(patient_metadata_nocytof_labels, 
                                                                                       align_widths = TRUE), 
                                                                         just = "left", 
                                                                         location = 0.25),
                                                   Sex = patient_metadata_nocytof$Sex,
                                                   Age = anno_barplot(patient_metadata_nocytof$Age, bar_width = 1, gp = gpar(fill = "#444444")),
                                                   BMI = anno_barplot(patient_metadata_nocytof$BMI, bar_width = 1, gp = gpar(fill = "#a69eb0")),
                                                   Group = patient_metadata_nocytof$Group_summarized, 
                                                   HIPEC = patient_metadata_nocytof$HIPEC,
                                                   `MS primary tumor` = patient_metadata_nocytof$MSprim,
                                                   `MS peritoneal metastasis` = patient_metadata_nocytof$MSPM,
                                                   # PM = patient_metadata_nocytof$PM,
                                                   `Pathological stage` = anno_simple(x = patient_metadata_nocytof$Cancer_stage, 
                                                                                      pch = patient_metadata_nocytof$Cancer_stage, 
                                                                                      col = c("I" = "#ffffff",
                                                                                              "II" = "#ffffff",
                                                                                              "III" = "#ffffff",
                                                                                              "IV" = "#ffffff"),
                                                                                      na_col = "#ffffff",
                                                                                      gp = gpar(col = "black")),
                                                   Omics = patient_metadata_nocytof$Omics,
                                                   annotation_name_side = "left", 
                                                   gp = gpar(col = "black"),
                                                   border = T,
                                                   na_col = "#ffffff",
                                                   col = list(Sex = c("Male" = "#3b55cd", 
                                                                      "Female" = "#cf3476"),
                                                              Group = c("HC" = "#91d1c2",
                                                                        "CRC-" = "#7f2284",
                                                                        "CRC+" = "#845422"),
                                                              Omics = c("CyTOF+scRNAseq" = "#009999",
                                                                        "scRNAseq" = "#FFFAA0",
                                                                        "CyTOF" = "#A7C7E7"),
                                                              `MS primary tumor` = c("MSI" = "#0000ff",
                                                                                     "MSS" = "#ff0000",
                                                                                     "Unknown" = "#a69eb0"),
                                                              `MS peritoneal metastasis` = c("MSI" = "#0000ff",
                                                                                             "MSS" = "#ff0000",
                                                                                             "Unknown" = "#a69eb0",
                                                                                             "No PM" = "#ffffff"),
                                                              HIPEC = c("Yes" = "#00BFC4",
                                                                        "No" = "#7CAE00")
                                                              # PM = c("No" = "#054f42",
                                                              #        "Yes" = "#AB2020")
                                                              )
                                                   )

pdf(width = 7.5, height = 7, file = heatmap_patient_metadata_nocytof_pdf)
plot(patient_metadata_nocytof_plot)
dev.off()

sessionInfo()