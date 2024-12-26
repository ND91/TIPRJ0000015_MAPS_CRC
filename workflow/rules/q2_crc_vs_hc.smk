############################
# Subsetting and combining #
############################

rule subsetting_hc_crcpmp_pbmc_pf:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crcpmp_pbmc_pf_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pbmc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pbmc_pf.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pbmc_pf_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_hc_crcpmp_pbmc_pf.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_crcpmp_pbmc_pf_seuratobject_rds}" &> "{log}"
    """

rule subsetting_hc_crcpmp_pf_allcells:
  input:
    curated_seuratobject_rds="output/curated/curated_SeuratObject.Rds",
  output:
    hc_crcpmp_pf_allcells_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_allcells_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf_allcells.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf_allcells_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_hc_crcpmp_pf_allcells.R "{input.curated_seuratobject_rds}" "{output.hc_crcpmp_pf_allcells_seuratobject_rds}" &> "{log}"
    """

rule subsetting_hc_crcpmp_pf:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crcpmp_pf_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_hc_crcpmp_pf.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_crcpmp_pf_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_pf_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    crcpmp_pf_macrophages_seuratobject_rds="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages_SeuratObject.Rds",
    crcpmp_pf_macrophages_cellmetadata_csv="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages/cellmetadata.csv",
    crcpmp_pf_macrophages_counts_mtx="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages/counts.mtx",
    crcpmp_pf_macrophages_features_csv="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages/features.csv",
    crcpmp_pf_macrophages_umap_csv="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages/umap.csv",
    crcpmp_pf_macrophages_pca_csv="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages/pca.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_crcpmp_pf_macrophages.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_crcpmp_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_crcpmp_pf_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.crcpmp_pf_macrophages_seuratobject_rds}" "{output.crcpmp_pf_macrophages_cellmetadata_csv}" "{output.crcpmp_pf_macrophages_counts_mtx}" "{output.crcpmp_pf_macrophages_features_csv}" "{output.crcpmp_pf_macrophages_umap_csv}" "{output.crcpmp_pf_macrophages_pca_csv}" &> "{log}"
    """

rule subsetting_hc_crcpmp_pbmc_pf_monocytes_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crcpmp_pbmc_pf_monocytes_macrophages_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pbmc_pf_monocytes_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pbmc_pf_monocytes_macrophages.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pbmc_pf_monocytes_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_hc_crcpmp_pbmc_pf_monocytes_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_crcpmp_pbmc_pf_monocytes_macrophages_seuratobject_rds}" &> "{log}"
    """

rule subsetting_hc_crcpmp_pf_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crcpmp_pf_macrophages_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages_SeuratObject.Rds",
    hc_crcpmp_pf_macrophages_cellmetadata_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/cellmetadata.csv",
    hc_crcpmp_pf_macrophages_counts_mtx="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/counts.mtx",
    hc_crcpmp_pf_macrophages_features_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/features.csv",
    hc_crcpmp_pf_macrophages_umap_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/umap.csv",
    hc_crcpmp_pf_macrophages_pca_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/pca.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf_macrophages.log",
  benchmark:
    "output/q2_crc_vs_hc/subsets/subsetting_hc_crcpmp_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q2_crc_vs_hc/subsetting_hc_crcpmp_pf_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_crcpmp_pf_macrophages_seuratobject_rds}" "{output.hc_crcpmp_pf_macrophages_cellmetadata_csv}" "{output.hc_crcpmp_pf_macrophages_counts_mtx}" "{output.hc_crcpmp_pf_macrophages_features_csv}" "{output.hc_crcpmp_pf_macrophages_umap_csv}" "{output.hc_crcpmp_pf_macrophages_pca_csv}" &> "{log}"
    """

############
# Analyses #
############

rule da_crcvhc:
  input:
    hc_crc_pbmc_pf_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pbmc_pf_SeuratObject.Rds",
    seuratDA_r="workflow/scripts/functions/seuratDA.R",
  output:
    dacs_list_rds="output/q2_crc_vs_hc/analyses/dacs_crcvhc_{tissue}_{comparison}_list.Rds",
    dacs_csv="output/q2_crc_vs_hc/analyses/dacs_crcvhc_{tissue}_{comparison}.csv",
    proportions_csv="output/q2_crc_vs_hc/analyses/proportions_crcvhc_{tissue}_{comparison}.csv",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/da_crcvhc_{tissue}_{comparison}.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/da_crcvhc_{tissue}_{comparison}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    tissue="{tissue}",
    comparison="{comparison}"
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/da_crcvhc.R "{input.hc_crc_pbmc_pf_seuratobject_rds}" "{input.seuratDA_r}" "{params.tissue}" "{params.comparison}" "{output.dacs_list_rds}" "{output.dacs_csv}" "{output.proportions_csv}" &> "{log}"
    """

rule de_crcvhc:
  input:
    hc_crc_pbmc_pf_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pbmc_pf_SeuratObject.Rds",
    seuratDE_r="workflow/scripts/functions/seuratDE.R",
  output:
    degs_list_rds="output/q2_crc_vs_hc/analyses/degs_crcvhc_{tissue}_{level}_list.Rds",
    degs_xlsx="output/q2_crc_vs_hc/analyses/degs_crcvhc_{tissue}_{level}_degs.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/de_crcvhc_{tissue}_{level}.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/de_crcvhc_{tissue}_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    tissue="{tissue}",
    level="{level}"
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/de_crcvhc.R "{input.hc_crc_pbmc_pf_seuratobject_rds}" "{input.seuratDE_r}" "{params.tissue}" "{params.level}" "{output.degs_list_rds}" "{output.degs_xlsx}" &> "{log}"
    """

rule de_fgsea_crcvhc:
  input:
    degs_list_rds="output/q2_crc_vs_hc/analyses/degs_crcvhc_{tissue}_{level}_list.Rds",
  output:
    fgsea_list_rds="output/q2_crc_vs_hc/analyses/fgsea_crcvhc_{tissue}_{level}_list.Rds",
    fgsea_pws_xlsx="output/q2_crc_vs_hc/analyses/fgsea_crcvhc_{tissue}_{level}_list.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/de_fgsea_crcvhc_{tissue}_{level}.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/de_fgsea_crcvhc_{tissue}_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/de_fgsea_crcvhc.R "{input.degs_list_rds}" "{output.fgsea_list_rds}" "{output.fgsea_pws_xlsx}" &> "{log}"
    """

rule de_fgsea_m1m2_azizi2018_crcvhc_pf_macrophages:
  input:
    degs_list_rds="output/q2_crc_vs_hc/analyses/degs_crcvhc_PF_manual_l4_list.Rds",
    azizi2018_xlsx="config/genes_of_interest/m1_m2_azizi2018.xlsx",
  output:
    fgsea_list_rds="output/q2_crc_vs_hc/analyses/fgsea_m1m2_azizi2018_crcvhc_pf_macrophages_list.Rds",
    fgsea_pws_xlsx="output/q2_crc_vs_hc/analyses/fgsea_m1m2_azizi2018_crcvhc_pf_macrophages_list.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/de_fgsea_m1m2_azizi2018_crcvhc_pf_macrophages.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/de_fgsea_m1m2_azizi2018_crcvhc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/de_fgsea_m1m2azizi2018_crcvhc.R "{input.degs_list_rds}" "{input.azizi2018_xlsx}" "{output.fgsea_list_rds}" "{output.fgsea_pws_xlsx}" &> "{log}"
    """

rule de_crc_pci_pf_macrophages:
  input:
    crc_pf_macrophages_seuratobject_rds="output/q2_crc_vs_hc/subsets/crcpmp_pf_macrophages_SeuratObject.Rds",
    seuratDE_r="workflow/scripts/functions/seuratDE.R",
  output:
    degs_list_rds="output/q2_crc_vs_hc/analyses/degs_crc_pci_pf_macrophages_{level}_list.Rds",
    degs_xlsx="output/q2_crc_vs_hc/analyses/degs_crc_pci_pf_macrophages_{level}_degs.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/de_crc_pci_pf_macrophages_{level}.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/de_crc_pci_pf_macrophages_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}"
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/de_crc_pci.R "{input.crc_pf_macrophages_seuratobject_rds}" "{input.seuratDE_r}" "{params.level}" "{output.degs_list_rds}" "{output.degs_xlsx}" &> "{log}"
    """

rule hc_crcpmp_pbmc_pf_monocytes_macrophages_trajectory_inference:
  input:
    hc_crcpmp_pbmc_pf_monocytes_macrophages_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pbmc_pf_monocytes_macrophages_SeuratObject.Rds",
  output:
    sce_ss_rds="output/q2_crc_vs_hc/analyses/hc_crcpmp_pbmc_pf_monocytes_macrophages_sce_ss.Rds",
  threads: 
    1
  conda:
    "../envs/r-trajectory.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/hc_crcpmp_pbmc_pf_monocytes_macrophages_trajectory_inference.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/hc_crcpmp_pbmc_pf_monocytes_macrophages_trajectory_inference_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/hc_crcpmp_pbmc_pf_monocytes_macrophages_trajectory_inference.R "{input.hc_crcpmp_pbmc_pf_monocytes_macrophages_seuratobject_rds}" "{output.sce_ss_rds}" &> "{log}"
    """
    
rule scvelo_analysis_hc_crcpmp_pf_macrophages:
  input:
    cellmetadata_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/cellmetadata.csv",
    counts_mtx="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/counts.mtx",
    features_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/features.csv",
    umap_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/umap.csv",
    pca_csv="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_macrophages/pca.csv",
    velocyto_curated_loom="output/curated/velocyto_curated.loom",
  output:
    scvelo_anndata_h5ad="output/q2_crc_vs_hc/analyses/hc_crcpmp_pf_macrophages_scvelo_anndata.h5ad",
    scvelo_cellmetadata_csv="output/q2_crc_vs_hc/analyses/hc_crcpmp_pf_macrophages_scvelo_cellmetadata.csv",
  threads: 
    1
  conda:
    "../envs/python-scvelo.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/scvelo_analysis_hc_crcpmp_pf_macrophages.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/scvelo_analysis_hc_crcpmp_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    python workflow/scripts/q1_pf_characterization/scvelo.py --cellmetadata_csv "{input.cellmetadata_csv}" --counts_mtx "{input.counts_mtx}" --features_csv "{input.features_csv}" --umap_csv "{input.umap_csv}" --pca_csv "{input.pca_csv}" --velocyto_curated_loom "{input.velocyto_curated_loom}" --scvelo_anndata_h5ad "{output.scvelo_anndata_h5ad}" --scvelo_anndata_cellmetadata_csv "{output.scvelo_cellmetadata_csv}" &> "{log}"
    """

rule hc_crcpmp_pf_nichenet:
  input:
    hc_crcpmp_pf_seuratobject_rds="output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_SeuratObject.Rds",
    degs_list_rds="output/q2_crc_vs_hc/analyses/degs_crcvhc_PF_manual_l3_list.Rds",
    lr_network_rds=config['lr_network'],
    ligand_target_matrix_rds=config['ligand_target_matrix'],
  output:
    nichenet_rds="output/q2_crc_vs_hc/analyses/hc_crcpmp_pf_nichenet_{celltype}.Rds",
  threads: 
    1
  conda:
    "../envs/r-nichenet.yaml",
  log:
    "output/q2_crc_vs_hc/analyses/hc_crcpmp_pf_nichenet_{celltype}.log",
  benchmark:
    "output/q2_crc_vs_hc/analyses/hc_crcpmp_pf_nichenet_{celltype}_benchmark.txt",
  resources:
    mem_mb=64000,
  params:
    celltype="{celltype}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/nichenet_crcvhc.R "{input.hc_crcpmp_pf_seuratobject_rds}" "{params.celltype}" "{input.lr_network_rds}" "{input.ligand_target_matrix_rds}" "{input.degs_list_rds}" "{output.nichenet_rds}" &> "{log}"
    """

###########
# Figures #
###########

rule fig_hc_crc_pmp_umap_pf_colcelltype_splitgroup:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_crc_pmp_umap_pf_colcelltype_splitgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q2_crc_vs_hc/fig_hc_crc_pmp_umap_pf_colcelltype_splitgroup.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_crc_pmp_umap_pf_colcelltype_splitgroup_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_umap_pf_colgroup:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_umap_pf_colgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_umap_pf_colgroup.pdf",
    hc_crc_pmp_umap_pf_colgroup_splitgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_umap_pf_colgroup_splitgroup.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/fig_hc_crc_pmp_{txdonor}_umap_pf_colgroup.log",
  benchmark:
    "output/figures/hc_crc/fig_hc_crc_pmp_{txdonor}_umap_pf_colgroup_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_umap_pf_colgroup.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{output.hc_crc_pmp_umap_pf_colgroup_pdf}" "{output.hc_crc_pmp_umap_pf_colgroup_splitgroup_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_boxplot_pbmc_pf_abundance:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_boxplot_pbmc_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_abundance_{comparison}.pdf",
    hc_crc_pmp_boxplot_pf_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_abundance_{comparison}.pdf",
    hc_crc_pmp_boxplot_pbmc_pf_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_boxplot_pbmc_pf_abundance.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.hc_crc_pmp_boxplot_pbmc_abundance_pdf}" "{output.hc_crc_pmp_boxplot_pf_abundance_pdf}" "{output.hc_crc_pmp_boxplot_pbmc_pf_abundance_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{output.hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_volcanoplot_pf_macrophage_l3_crcpmpvhc:
  input:
    deseq2_list_rds="output/analyses/hc_crc_pmp_{txdonor}_pf_macrophages_l3_crcpmpvhc_deseq2_list.Rds",
  output:
    volcanoplot_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_crcpmpvhc.pdf",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_macrophage_l3_crcpmpvhc.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_macrophage_l3_crcpmpvhc_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_volcanoplot_pf_macrophages_l3_crcpmpvhc_de.R "{input.deseq2_list_rds}" "{output.volcanoplot_pdf}" &> "{log}"
    """
