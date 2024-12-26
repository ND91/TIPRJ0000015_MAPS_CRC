# Part 1: PF of healthy individuals

## Subsetting and combining cells

rule subsetting_hc_pbmc_pf:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pbmc_pf.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pbmc_pf_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pbmc_pf.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_pbmc_pf_seuratobject_rds}" &> "{log}"
    """
    
rule subsetting_hc_pbmc:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_pbmc_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pbmc.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pbmc_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pbmc.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_pbmc_seuratobject_rds}" &> "{log}"
    """
    
rule subsetting_hc_pf:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_pf_seuratobject_rds}" &> "{log}"
    """

rule subsetting_hc_pf_subset:
  input:
    hc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds",
  output:
    hc_pf_subset_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_{lineage}_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_{lineage}.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_{lineage}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    lineage="{lineage}",
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf_subset.R "{input.hc_pf_seuratobject_rds}" "{params.lineage}" "{output.hc_pf_subset_seuratobject_rds}" &> "{log}"
    """

rule subsetting_hc_pf_macrophages:
  input:
    hc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds",
  output:
    hc_pf_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_SeuratObject.Rds",
    hc_pf_macrophages_cellmetadata_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/cellmetadata.csv",
    hc_pf_macrophages_counts_mtx="output/q1_pf_characterization/subsets/hc_pf_macrophages/counts.mtx",
    hc_pf_macrophages_features_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/features.csv",
    hc_pf_macrophages_umap_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/umap.csv",
    hc_pf_macrophages_pca_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/pca.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf_macrophages.R "{input.hc_pf_seuratobject_rds}" "{output.hc_pf_macrophages_seuratobject_rds}" "{output.hc_pf_macrophages_cellmetadata_csv}" "{output.hc_pf_macrophages_counts_mtx}" "{output.hc_pf_macrophages_features_csv}" "{output.hc_pf_macrophages_umap_csv}" "{output.hc_pf_macrophages_pca_csv}" &> "{log}"
    """

rule combining_pbmc_pf_liver_colon:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    liver_seuratobject_rds="resources/liver/liver_reannotated_seuratobject.rds",
    colon_seuratobject_rds="resources/colon/colon_reannotated_seuratobject.rds",
  output:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/combining_pbmc_pf_liver_colon.log",
  benchmark:
    "output/q1_pf_characterization/subsets/combining_pbmc_pf_liver_colon_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/combine_pbmc_pf_liver_colon.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.liver_seuratobject_rds}" "{input.colon_seuratobject_rds}" "{output.pbmc_pf_liver_colon_seuratobject_rds}" &> "{log}"
    """   

rule subsetting_liver_mnp:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    liver_mnp_seuratobject_rds="output/q1_pf_characterization/subsets/liver_mnp_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_liver_mnp.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_liver_mnp_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_liver_mnp.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.liver_mnp_seuratobject_rds}" &> "{log}"
    """

rule subsetting_liver_macrophages:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    liver_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/liver_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_liver_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_liver_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_liver_macrophages.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.liver_macrophages_seuratobject_rds}" &> "{log}"
    """

rule subsetting_colon_mnp:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    colon_mnp_seuratobject_rds="output/q1_pf_characterization/subsets/colon_mnp_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_colon_mnp.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_colon_mnp_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_colon_mnp.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.colon_mnp_seuratobject_rds}" &> "{log}"
    """

rule subsetting_colon_macrophages:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    colon_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/colon_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_colon_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_colon_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_colon_macrophages.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.colon_macrophages_seuratobject_rds}" &> "{log}"
    """

## Analyses

rule markergenes_hc_pf_macrophages:
  input:
    hc_pf_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_SeuratObject.Rds",
  output:
    hc_pf_macrophages_marker_list_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_markergenes_list.Rds",
  threads: 
    8
  conda:
    "../envs/r-ucell.yaml",
  log:
    "output/q1_pf_characterization/analyses/markergenes_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/markergenes_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/markergenes_hc_pf_macrophages_manual_l4.R "{input.hc_pf_macrophages_seuratobject_rds}" "{output.hc_pf_macrophages_marker_list_rds}" &> "{log}"
    """

rule markerproteins_hc_pf_manual_l2l3:
  input:
    hc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds",
  output:
    hc_pf_marker_list_rds="output/q1_pf_characterization/analyses/hc_pf_manual_l2l3_markerproteins_list.Rds",
  threads: 
    8
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/analyses/markerproteins_hc_pf_manual_l2l3.log",
  benchmark:
    "output/q1_pf_characterization/analyses/markerproteins_hc_pf_manual_l2l3_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/markerproteins_hc_pf_manual_l2l3.R "{input.hc_pf_seuratobject_rds}" "{output.hc_pf_marker_list_rds}" &> "{log}"
    """

rule markerproteins_hc_pf_macrophages:
  input:
    hc_pf_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_SeuratObject.Rds",
  output:
    hc_pf_macrophages_marker_list_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_markerproteins_list.Rds",
  threads: 
    8
  conda:
    "../envs/r-ucell.yaml",
  log:
    "output/q1_pf_characterization/analyses/markerproteins_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/markerproteins_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/markerproteins_hc_pf_macrophages_manual_l4.R "{input.hc_pf_macrophages_seuratobject_rds}" "{output.hc_pf_macrophages_marker_list_rds}" &> "{log}"
    """

rule markergenes_fgsea_hc_pf_macrophages:
  input:
    hc_pf_macrophages_marker_list_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_marker_list.Rds",
  output:
    fgsea_list_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_marker_fgsea_list.Rds",
    fgsea_pws_xlsx="output/q1_pf_characterization/analyses/hc_pf_macrophages_marker_fgsea_list.xlsx",
  threads: 
    8
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/markergenes_fgsea_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/markergenes_fgsea_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/markergenes_fgsea_hc_pf_macrophages_manual_l4.R "{input.hc_pf_macrophages_marker_list_rds}" "{output.fgsea_list_rds}" "{output.fgsea_pws_xlsx}" &> "{log}"
    """

rule da_hc_pfvpbmc:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    seuratDA_r="workflow/scripts/functions/seuratDA.R",
  output:
    dacs_list_rds="output/q1_pf_characterization/analyses/hc_pfvpbmc_{comparison}_dacs_list.Rds",
    dacs_csv="output/q1_pf_characterization/analyses/hc_pfvpbmc_{comparison}_dacs.csv",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/da_hc_pfvpbmc_{comparison}.log",
  benchmark:
    "output/q1_pf_characterization/analyses/da_hc_pfvpbmc_{comparison}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/da_hc_pfvpbmc.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.seuratDA_r}" "{params.comparison}" "{output.dacs_list_rds}" "{output.dacs_csv}" &> "{log}"
    """

rule de_hc_pfvpbmc:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    seuratDE_r="workflow/scripts/functions/seuratDE.R",
  output:
    deseq2_list_rds="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_deseq2_list.Rds",
    degs_xlsx="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_degs.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/de_hc_pfvpbmc_{level}.log",
  benchmark:
    "output/q1_pf_characterization/analyses/de_hc_pfvpbmc_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/de_hc_pfvpbmc.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.seuratDE_r}" "{params.level}" "{output.deseq2_list_rds}" "{output.degs_xlsx}" &> "{log}"
    """
    
rule de_fgsea_hc_pfvpbmc:
  input:
    deseq2_list_rds="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_deseq2_list.Rds",
  output:
    fgsea_list_rds="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_fgsea_list.Rds",
    fgsea_pws_xlsx="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_fgsea_list.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/de_fgsea_hc_pfvpbmc_{level}.log",
  benchmark:
    "output/q1_pf_characterization/analyses/de_fgsea_hc_pfvpbmc_{level}_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/de_fgsea_hc_pfvpbmc.R "{input.deseq2_list_rds}" "{output.fgsea_list_rds}" "{output.fgsea_pws_xlsx}" &> "{log}"
    """

rule tam_classification_hc_pf_macrophages:
  input:
    hc_pf_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_SeuratObject.Rds",
    tam_markers_xlsx=config['tam_markers'],
  output:
    hc_pf_macrophages_tamannotation_scores_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_tamannotation_scores.Rds",
    hc_pf_macrophages_tamannotation_ranks_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_tamannotation_ranks.Rds",
  threads: 
    8
  conda:
    "../envs/r-ucell.yaml",
  log:
    "output/q1_pf_characterization/analyses/tam_classification_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/tam_classification_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/tam_classification_hc_pf_macrophages.R "{input.hc_pf_macrophages_seuratobject_rds}" "{input.tam_markers_xlsx}" "{threads}" "{output.hc_pf_macrophages_tamannotation_scores_rds}" "{output.hc_pf_macrophages_tamannotation_ranks_rds}" &> "{log}"
    """

# Part 2: PF of PM-CRC patients

## Subsetting and combining cells

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

## Analyses

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

# Part 3: TX of PM-CRC patients

## Subsetting and combining cells 

rule subsetting_hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_tx_immune:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
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
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_tx_immune.R "{input.seurat_curated_rds}" "{output.crcpmp_tx_immune_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_tx_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    crcpmp_tx_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_tx_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.crcpmp_tx_macrophages_seuratobject_rds}" &> "{log}"
    """

rule subsetting_crcpmp_tx_monocytes_macrophages:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    crcpmp_tx_monocytes_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages_SeuratObject.Rds",
    crcpmp_tx_monocytes_macrophages_cellmetadata_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/cellmetadata.csv",
    crcpmp_tx_monocytes_macrophages_counts_mtx="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/counts.mtx",
    crcpmp_tx_monocytes_macrophages_features_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/features.csv",
    crcpmp_tx_monocytes_macrophages_umap_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/umap.csv",
    crcpmp_tx_monocytes_macrophages_pca_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/pca.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_monocytes_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/subsets/subsetting_crcpmp_tx_monocytes_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/subsetting_crcpmp_tx_monocytes_macrophages.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.crcpmp_tx_monocytes_macrophages_seuratobject_rds}" "{output.crcpmp_tx_monocytes_macrophages_cellmetadata_csv}" "{output.crcpmp_tx_monocytes_macrophages_counts_mtx}" "{output.crcpmp_tx_monocytes_macrophages_features_csv}" "{output.crcpmp_tx_monocytes_macrophages_umap_csv}" "{output.crcpmp_tx_monocytes_macrophages_pca_csv}" &> "{log}"
    """

## Analyses

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

rule tam_classification_crcpmp_pf_tx_macrophages:
  input:
    crcpmp_pf_tx_paired_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_SeuratObject.Rds",
    tam_markers_xlsx=config['tam_markers'],
  output:
    crcpmp_pf_tx_paired_macrophages_tamannotation_seuratobject_rds="output/q3_pm_tx_characterization/analyses/crcpmp_pf_tx_paired_macrophages_tamannotation.Rds",
  threads: 
    8
  conda:
    "../envs/r-ucell.yaml",
  log:
    "output/q3_pm_tx_characterization/analyses/tam_classification_crcpmp_pf_tx_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/analyses/tam_classification_crcpmp_pf_tx_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q3_pm_tx_characterization/tam_classification_crcpmp_pf_tx_macrophages.R "{input.crcpmp_pf_tx_paired_macrophages_seuratobject_rds}" "{input.tam_markers_xlsx}" "{threads}" "{output.crcpmp_pf_tx_paired_macrophages_tamannotation_seuratobject_rds}" &> "{log}"
    """
    
rule scvelo_analysis_crcpmp_tx_monocytes_macrophages:
  input:
    cellmetadata_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/cellmetadata.csv",
    counts_mtx="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/counts.mtx",
    features_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/features.csv",
    umap_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/umap.csv",
    pca_csv="output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages/pca.csv",
    velocyto_curated_loom="output/curated/velocyto_curated.loom",
  output:
    scvelo_anndata_h5ad="output/q3_pm_tx_characterization/analyses/crcpmp_tx_monocytes_macrophages_scvelo_anndata.h5ad",
    scvelo_cellmetadata_csv="output/q3_pm_tx_characterization/analyses/crcpmp_tx_monocytes_macrophages_scvelo_cellmetadata.csv",
  threads: 
    1
  conda:
    "../envs/python-scvelo.yaml",
  log:
    "output/q3_pm_tx_characterization/analyses/scvelo_analysis_crcpmp_tx_monocytes_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/analyses/scvelo_analysis_crcpmp_tx_monocytes_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    python workflow/scripts/q3_pm_tx_characterization/scvelo_crcpmp_tx_monocytes_macrophages.py --cellmetadata_csv "{input.cellmetadata_csv}" --counts_mtx "{input.counts_mtx}" --features_csv "{input.features_csv}" --umap_csv "{input.umap_csv}" --pca_csv "{input.pca_csv}" --velocyto_curated_loom "{input.velocyto_curated_loom}" --scvelo_anndata_h5ad "{output.scvelo_anndata_h5ad}" --scvelo_anndata_cellmetadata_csv "{output.scvelo_cellmetadata_csv}" &> "{log}"
    """

rule trajectory_analysis_crcpmp_tx_macrophages:
  input:
    crcpmp_tx_macrophages_seuratobject_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_SeuratObject.Rds",
  output:
    crcpmp_tx_macrophages_trajectory_sce_rds="output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_trajectory_sce.Rds",
    aggregated_lines_csv="output/q3_pm_tx_characterization/analyses/crcpmp_tx_macrophages_trajectory_aggregated_lines.csv",
    tscan_rootmacrophagesvcan_rds="output/q3_pm_tx_characterization/analyses/crcpmp_tx_macrophages_trajectory_tscan_rootmacrophagesvcan.Rds",
  threads: 
    1
  conda:
    "../envs/r-trajectory.yaml",
  log:
    "output/q3_pm_tx_characterization/analyses/trajectory_analysis_crcpmp_tx_macrophages.log",
  benchmark:
    "output/q3_pm_tx_characterization/analyses/trajectory_analysis_crcpmp_tx_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q3_pm_tx_characterization/trajectory_analysis_crcpmp_tx_macrophages.R "{input.crcpmp_tx_macrophages_seuratobject_rds}" "{output.crcpmp_tx_macrophages_trajectory_sce_rds}" "{output.aggregated_lines_csv}" "{output.tscan_rootmacrophagesvcan_rds}" &> "{log}"
    """
