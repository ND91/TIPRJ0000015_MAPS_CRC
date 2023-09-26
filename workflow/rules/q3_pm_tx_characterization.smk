##############
# Subsetting #
##############

rule subsetting_crcpmp_pbmc_pf_tx_paired:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    crcpmp_pbmc_pf_tx_paired_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pbmc_pf_tx_paired_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pbmc_pf_tx_paired.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pbmc_pf_tx_paired_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_pbmc_pf_tx_paired.R "{input.seurat_curated_rds}" "{output.crcpmp_pbmc_pf_tx_paired_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_pbmc_pf_tx_paired_anndata:
  input:
    curated_loom="output/curated/velocyto_curated.loom",
  output:
    crcpmp_pbmc_pf_tx_paired_anndata_h5ad="output/q3_pm_tx_characterization/subsets/crcpmp_pbmc_pf_tx_paired_anndata.h5ad",
  threads: 
    1
  conda:
    "../envs/scvelo.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pbmc_pf_tx_paired_anndata.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pbmc_pf_tx_paired_anndata_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_pbmc_pf_tx_paired_anndata.py --curated_loom "{input.curated_loom}" --subset_adata_h5ad "{output.crcpmp_pbmc_pf_tx_paired_anndata_h5ad}" &> "{log}"
    """

rule subsetting_crcpmp_pf_tx_paired:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    crc_pf_tx_paired_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crc_pf_tx_paired_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pf_tx_paired.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pf_tx_paired_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_pf_tx_paired.R "{input.seurat_curated_rds}" "{output.crc_pf_tx_paired_seuratobject_rds}" &> "{log}"
    """
    
rule subsetting_crcpmp_pf_tx_paired_immune:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    crcpmp_pf_tx_paired_immune_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_immune_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pf_tx_paired_immune.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_pf_tx_paired_immune_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_pf_tx_paired_immune.R "{input.seurat_curated_rds}" "{output.crcpmp_pf_tx_paired_immune_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_tx_immune:
  input:
    crcpmp_pf_tx_paired_immune_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_immune_SeuratObject.Rds",
  output:
    crcpmp_tx_immune_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_immune_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_immune.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_immune_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_tx.R "{input.crcpmp_pf_tx_paired_immune_seuratobject_rds}" "{output.crcpmp_tx_immune_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_pbmc_pf_tx_t_paired_cdr3:
  input:
    airr_merged_csv="output/trust4/trust4_barcode_annotated_airr.csv",
  output:
    crcpmp_pbmc_pf_tx_t_paired_cdr3_airr_csv="output/q3_pm_tx_characterization/trust4/crcpmp_pbmc_pf_tx_t_paired_airr.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/trust4/subsetting_crcpmp_pbmc_pf_tx_t_paired_cdr3.log",
  benchmark:
    "output/q3_pm_tx_characterization/trust4/subsetting_crcpmp_pbmc_pf_tx_t_paired_cdr3_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/subsetting_crcpmp_pbmc_pf_tx_t_paired_cdr3_subsetting.R "{input.airr_merged_csv}" "{output.crcpmp_pbmc_pf_tx_t_paired_cdr3_airr_csv}" &> "{log}"
    """
    
rule subsetting_crcpmp_pf_tx_paired_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    crcpmp_pf_tx_paired_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_hc_crcpmp_pf_tx_paired_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_hc_crcpmp_pf_tx_paired_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_pf_tx_paired_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.crcpmp_pf_tx_paired_macrophages_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_tx_macrophages_subsetting:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    crcpmp_tx_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_subsetting.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_tx_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.crcpmp_tx_macrophages_seuratobject_rds}" &> "{log}"
    """

############
# Analyses #
############

rule trajectory_analysis_crcpmp_pf_tx_paired_macrophages:
  input:
    crcpmp_pf_tx_paired_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_SeuratObject.Rds",
  output:
    crcpmp_pf_tx_paired_macrophages_trajectory_sce_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_trajectory_sce.Rds",
    aggregated_lines_csv="output/q3_pm_tx_characterization/analyses/crcpmp_pf_tx_paired_macrophages_trajectory_aggregated_lines.csv",
    tscan_rootcl1_rds="output/q3_pm_tx_characterization/analyses/crcpmp_pf_tx_paired_macrophages_trajectory_tscan_rootcl1.Rds",
  threads: 
    1
  conda:
    "../envs/r-trajectory.yaml",
  log:
    "output/q3_pm_tx_characterization/analyses/trajectory_analysis_crcpmp_pf_tx_paired_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/analyses/trajectory_analysis_crcpmp_pf_tx_paired_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q3_pm_tx_characterization/trajectory_analysis_crcpmp_pf_tx_paired_macrophages.R "{input.crcpmp_pf_tx_paired_macrophages_seuratobject_rds}" "{output.crcpmp_pf_tx_paired_macrophages_trajectory_sce_rds}" "{output.aggregated_lines_csv}" "{output.tscan_rootcl1_rds}" &> "{log}"
    """

# rule hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_trajectory_inference:
#   input:
#     hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_SeuratObject.Rds",
#   output:
#     sce_ss_rds="output/q3_pm_tx_characterization/analyses/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_sce_ss.Rds",
#   threads: 
#     1
#   conda:
#     "../envs/r-trajectory.yaml",
#   log:
#     "output/q3_pm_tx_characterization/analyses/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_trajectory_inference.log",
#   benchmark:
#     "output/q3_pm_tx_characterization/analyses/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_trajectory_inference_benchmark.txt",
#   resources:
#     mem_mb=64000,
#   shell:
#     """
#     Rscript --vanilla workflow/scripts/q3_pm_tx_characterization/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_trajectory_inference.R "{input.hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_seuratobject_rds}" "{output.sce_ss_rds}" &> "{log}"
#     """
    
rule de_crcpmp_txvpf:
  input:
    crc_pf_tx_paired_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crc_pf_tx_paired_SeuratObject.Rds",
    seuratDE_r="workflow/scripts/functions/seuratDE.R",
  output:
    deseq2_list_rds="output/q3_pm_tx_characterization/analyses/crcpmp_txvpf_{level}_deseq2_list.Rds",
    degs_xlsx="output/q3_pm_tx_characterization/analyses/crcpmp_txvpf_{level}_degs.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q3_pm_tx_characterization/analyses/de_crcpmp_txvpf_{level}.log",
  benchmark:
    "output/q3_pm_tx_characterization/analyses/de_crcpmp_txvpf_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q3_pm_tx_characterization/de_crcpmp_txvpf.R "{input.crc_pf_tx_paired_seuratobject_rds}" "{input.seuratDE_r}" "{params.level}" "{output.deseq2_list_rds}" "{output.degs_xlsx}" &> "{log}"
    """

###########
# Figures #
###########

rule fig_crcpmp_patient_metadata:
  input:
    patient_metadata_xlsx=config["patient_metadata"],
  output:
    crc_heatmap_patient_metadata_pdf="output/figures/crc/crc_heatmap_patient_metadata.pdf",
    crc_heatmap_patient_metadata_nocytof_pdf="output/figures/crc/crc_heatmap_patient_metadata_nocytof.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/crc/crc_patient_metadata.log",
  benchmark:
    "output/figures/crc/crc_patient_metadata_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_heatmap_patient_metadata.R "{input.patient_metadata_xlsx}" "{output.crc_heatmap_patient_metadata_pdf}" "{output.crc_heatmap_patient_metadata_nocytof_pdf}" &> "{log}"
    """ 

rule fig_crcpmp_umap_pbmc_pf_tx_colcelltype_splittissue:
  input:
    crc_pbmc_pf_tx_paired_seuratobject_rds="output/subsets/crc_pbmc_pf_tx_paired_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_umap_pbmc_pf_tx_paired_colcelltype_splittissue_pdf="output/figures/crc/crc_umap_pbmc_pf_tx_paired_col{level}_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_umap_pbmc_pf_tx_col{level}_splittissue.log",
  benchmark:
    "output/figures/crc/crc_umap_pbmc_pf_tx_col{level}_splittissue_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_umap_pbmc_pf_tx_colcelltype_splittissue.R "{input.crc_pbmc_pf_tx_paired_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_umap_pbmc_pf_tx_paired_colcelltype_splittissue_pdf}" &> "{log}"
    """

rule fig_crcpmp_boxplot_tx_abundance:
  input:
    crcpmp_tx_seuratobject_rds="output/subsets/crcpmp_tx_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_boxplot_tx_abundance_top5_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_top5.pdf",
    crc_boxplot_tx_abundance_top10_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_top10.pdf",
    crc_boxplot_tx_abundance_all_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_all.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_tx_abundance_{level}.log",
  benchmark:
    "output/figures/crc/crc_boxplot_tx_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_boxplot_tx_abundance.R "{input.crc_tx_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_boxplot_tx_abundance_top5_pdf}" "{output.crc_boxplot_tx_abundance_top10_pdf}" "{output.crc_boxplot_tx_abundance_all_pdf}" &> "{log}"
    """

rule fig_crcpmp_boxplot_myeloid_tx_abundance:
  input:
    crc_myeloid_tx_seuratobject_rds="output/subsets/crc_tx_myeloid_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_boxplot_myeloid_tx_abundance_top5_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_top5.pdf",
    crc_boxplot_myeloid_tx_abundance_top10_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_top10.pdf",
    crc_boxplot_myeloid_tx_abundance_all_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_all.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}.log",
  benchmark:
    "output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_boxplot_myeloid_tx_abundance.R "{input.crc_myeloid_tx_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_boxplot_myeloid_tx_abundance_top5_pdf}" "{output.crc_boxplot_myeloid_tx_abundance_top10_pdf}" "{output.crc_boxplot_myeloid_tx_abundance_all_pdf}" &> "{log}"
    """

rule fig_crcpmp_scatterplot_pf_tx_macrophage_l4rl2_pci:
  input:
    crc_pbmc_pf_tx_paired_seuratobject_rds="output/subsets/crc_pbmc_pf_tx_paired_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_scatterplot_pf_tx_macrophage_l4rl2_pci_pdf="output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci.log",
  benchmark:
    "output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_scatterplot_macrophage_l4_pci.R "{input.crc_pbmc_pf_tx_paired_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{output.crc_scatterplot_pf_tx_macrophage_l4rl2_pci_pdf}" &> "{log}"
    """

rule fig_crcpmp_t_clonotype_abundance:
  input:
    crc_pbmc_pf_tx_t_paired_cdr3_airr_csv="output/trust4/crc_pbmc_pf_tx_t_paired_airr.csv",
  output:
    crc_upset_t_clonotype_abundance_pdf="output/figures/crc/crc_upset_t_clonotype_abundance_{patient}.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/crc/crc_t_clonotype_abundance_{patient}.log",
  benchmark:
    "output/figures/crc/crc_t_clonotype_abundance_{patient}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    patient="{patient}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_t_clonotype_abundance.R "{input.crc_pbmc_pf_tx_t_paired_cdr3_airr_csv}" "{params.patient}" "{output.crc_upset_t_clonotype_abundance_pdf}" &> "{log}"
    """

rule fig_crcpmp_boxplot_t_pairwise_clonotype_abundance:
  input:
    crc_pbmc_pf_tx_t_paired_cdr3_airr_csv="output/trust4/crc_pbmc_pf_tx_t_paired_airr.csv",
  output:
    crc_upset_t_clonotype_abundance_pdf="output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance.pdf",
    crc_upset_t_clonotype_abundance_patanno_pdf="output/figures/crc/crc_boxplot_t_pairwise_clonotype_patanno_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance.log",
  benchmark:
    "output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crcpmp_boxplot_t_pairwise_clonotype_abundance.R "{input.crc_pbmc_pf_tx_t_paired_cdr3_airr_csv}" "{output.crc_upset_t_clonotype_abundance_pdf}" "{output.crc_upset_t_clonotype_abundance_patanno_pdf}" &> "{log}"
    """










rule fig_boxplot_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/subsets/tx_SeuratObject.Rds",
  output:
    boxplot_tx_l3rl2_df_csv="output/figures/fig_boxplot_tx_l3rl2.csv",
    boxplot_tx_l3rl2_svg="output/figures/fig_boxplot_tx_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_tx_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_tx_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_tx_l3rl2.R {input.tx_seuratobject_rds} {output.boxplot_tx_l3rl2_df_csv} {output.boxplot_tx_l3rl2_svg} &> {log}
    """

rule fig_tsne_pbmc_pf_tx_tsne:
  input:
    pbmc_pf_tx_paired_seuratobject_rds="output/subsets/pbmc_pf_tx_paired_SeuratObject.Rds",
  output:
    tsne_pbmc_pf_tx_df_csv="output/figures/tsne_pbmc_pf_tx.csv",
    tsne_pbmc_pf_tx_fig1_svg="output/figures/fig_tsne_pbmc_pf_tx_fig1.svg",
    tsne_pbmc_pf_tx_fig2_svg="output/figures/fig_tsne_pbmc_pf_tx_fig2.svg",
    tsne_pbmc_pf_tx_fig3_svg="output/figures/fig_tsne_pbmc_pf_tx_fig3.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_pbmc_pf_tx_tsne.log",
  benchmark:
    "output/figures/fig_pbmc_pf_tx_tsne_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_tsne_pbmc_pf_tx.R {input.pbmc_pf_tx_paired_seuratobject_rds} {output.tsne_pbmc_pf_tx_df_csv} {output.tsne_pbmc_pf_tx_fig1_svg} {output.tsne_pbmc_pf_tx_fig2_svg} {output.tsne_pbmc_pf_tx_fig3_svg} &> {log}
    """

rule fig_boxplot_pf_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/subsets/tx_SeuratObject.Rds",
  output:
    boxplot_pf_tx_l3rl2_df_csv="output/figures/fig_boxplot_pf_tx_l3rl2.csv",
    boxplot_pf_tx_l3rl2_svg="output/figures/fig_boxplot_pf_tx_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_pf_tx_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_pf_tx_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_pf_tx_l3rl2.R "{input.tx_seuratobject_rds}" "{output.boxplot_pf_tx_l3rl2_df_csv}" "{output.boxplot_pf_tx_l3rl2_svg}" &> "{log}"
    """

rule fig_boxplot_pbmc_pf_l3rl2:
  input:
    pbmc_pf_tx_seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
  output:
    boxplot_pbmc_pf_l3rl2_df_csv="output/figures/fig_boxplot_pbmc_pf_l3rl2.csv",
    boxplot_pbmc_pf_l3rl2_svg="output/figures/fig_boxplot_pbmc_pf_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_pbmc_pf_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_pbmc_pf_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_pbmc_pf_l3rl2.R "{input.pbmc_pf_tx_seuratobject_rds}" "{output.boxplot_pbmc_pf_l3rl2_df_csv}" "{output.boxplot_pbmc_pf_l3rl2_svg}" &> "{log}"
    """

rule fig_heatmap_cd4t_markers:
  input:
    pbmc_pf_tx_seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    markers_csv=config["markers"],
  output:
    heatmap_cd4t_markers_svg="output/figures/fig_heatmap_cd4t_markers.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_heatmap_cd4t_markers.log",
  benchmark:
    "output/figures/fig_heatmap_cd4t_markers_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_heatmap_markers_cd4t.R {input.pbmc_pf_tx_seuratobject_rds} {input.markers_csv} {output.heatmap_cd4t_markers_svg} &> {log}
    """
