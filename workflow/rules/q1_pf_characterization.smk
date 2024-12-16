############################
# Subsetting and combining #
############################

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
    hc_pf_subset_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_{subset}_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_{subset}.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_hc_pf_{subset}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    subset="{subset}",
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf_subset.R "{input.hc_pf_seuratobject_rds}" "{params.subset}" "{output.hc_pf_subset_seuratobject_rds}" &> "{log}"
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

rule subsetting_pf_liver_colon_mnp:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    pbmc_pf_liver_colon_mnp_seuratobject_rds="output/q1_pf_characterization/subsets/pf_liver_colon_mnp_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_pf_liver_colon_mnp.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_pf_liver_colon_mnp_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf_liver_colon_mnp.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.pbmc_pf_liver_colon_mnp_seuratobject_rds}" &> "{log}"
    """

rule subsetting_pf_liver_colon_macrophages:
  input:
    pbmc_pf_liver_colon_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
  output:
    pbmc_pf_liver_colon_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/pf_liver_colon_macrophages_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/subsets/subsetting_pf_liver_colon_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/subsets/subsetting_pf_liver_colon_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/q1_pf_characterization/subsetting_hc_pf_liver_colon_macrophages.R "{input.pbmc_pf_liver_colon_seuratobject_rds}" "{output.pbmc_pf_liver_colon_macrophages_seuratobject_rds}" &> "{log}"
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

############
# Analyses #
############

rule trajectory_analysis_hc_pf_macrophages:
  input:
    hc_pf_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_SeuratObject.Rds",
  output:
    hc_pf_macrophages_trajectory_sce_rds="output/q1_pf_characterization/subsets/hc_pf_macrophages_trajectory_sce.Rds",
    aggregated_lines_csv="output/q1_pf_characterization/analyses/hc_pf_macrophages_trajectory_aggregated_lines.csv",
    tscan_rootmacrophagesvcan_rds="output/q1_pf_characterization/analyses/hc_pf_macrophages_trajectory_tscan_rootmacrophagesvcan.Rds",
  threads: 
    1
  conda:
    "../envs/r-trajectory.yaml",
  log:
    "output/q1_pf_characterization/analyses/trajectory_analysis_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/trajectory_analysis_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/trajectory_analysis_pf_macrophages.R "{input.hc_pf_macrophages_seuratobject_rds}" "{output.hc_pf_macrophages_trajectory_sce_rds}" "{output.aggregated_lines_csv}" "{output.tscan_rootmacrophagesvcan_rds}" &> "{log}"
    """

rule scvelo_analysis_hc_pf_macrophages:
  input:
    cellmetadata_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/cellmetadata.csv",
    counts_mtx="output/q1_pf_characterization/subsets/hc_pf_macrophages/counts.mtx",
    features_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/features.csv",
    umap_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/umap.csv",
    pca_csv="output/q1_pf_characterization/subsets/hc_pf_macrophages/pca.csv",
    velocyto_curated_loom="output/curated/velocyto_curated.loom",
  output:
    scvelo_anndata_h5ad="output/q1_pf_characterization/analyses/hc_pf_macrophages_scvelo_anndata.h5ad",
    scvelo_cellmetadata_csv="output/q1_pf_characterization/analyses/hc_pf_macrophages_scvelo_cellmetadata.csv",
  threads: 
    1
  conda:
    "../envs/python-scvelo.yaml",
  log:
    "output/q1_pf_characterization/analyses/scvelo_analysis_hc_pf_macrophages.log",
  benchmark:
    "output/q1_pf_characterization/analyses/scvelo_analysis_hc_pf_macrophages_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    python workflow/scripts/q1_pf_characterization/scvelo.py --cellmetadata_csv "{input.cellmetadata_csv}" --counts_mtx "{input.counts_mtx}" --features_csv "{input.features_csv}" --umap_csv "{input.umap_csv}" --pca_csv "{input.pca_csv}" --velocyto_curated_loom "{input.velocyto_curated_loom}" --scvelo_anndata_h5ad "{output.scvelo_anndata_h5ad}" --scvelo_anndata_cellmetadata_csv "{output.scvelo_cellmetadata_csv}" &> "{log}"
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

rule da_hc_pfvlivervcolon:
  input:
    pbmc_pf_liver_colon_subset_seuratobject_rds="output/q1_pf_characterization/subsets/pbmc_pf_liver_colon_SeuratObject.Rds",
    seuratDA_r="workflow/scripts/functions/seuratDA.R",
  output:
    dacs_anova_rds="output/q1_pf_characterization/analyses/hc_pfvlivervcolon_{comparison}_dacs_list.Rds",
    dacs_anova_csv="output/q1_pf_characterization/analyses/hc_pfvlivervcolon_{comparison}_dacs.csv",
    dacs_livervpf_rds="output/q1_pf_characterization/analyses/hc_livervpf_{comparison}_dacs_list.Rds",
    dacs_livervpf_csv="output/q1_pf_characterization/analyses/hc_livervpf_{comparison}_dacs.csv",
    dacs_colonvpf_rds="output/q1_pf_characterization/analyses/hc_colonvpf_{comparison}_dacs_list.Rds",
    dacs_colonvpf_csv="output/q1_pf_characterization/analyses/hc_colonvpf_{comparison}_dacs.csv",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/da_hc_pfvlivervcolon_{comparison}.log",
  benchmark:
    "output/q1_pf_characterization/analyses/da_hc_pfvlivervcolon_{comparison}_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/da_hc_pfvlivervcolon.R "{input.pbmc_pf_liver_colon_subset_seuratobject_rds}" "{input.seuratDA_r}" "{params.comparison}" "{output.dacs_anova_rds}" "{output.dacs_anova_csv}" "{output.dacs_livervpf_rds}" "{output.dacs_livervpf_csv}" "{output.dacs_colonvpf_rds}" "{output.dacs_colonvpf_csv}" &> "{log}"
    """

rule de_hc_pfvlivervcolon:
  input:
    pbmc_pf_liver_colon_macrophages_seuratobject_rds="output/q1_pf_characterization/subsets/pf_liver_colon_macrophages_SeuratObject.Rds",
    seuratDE_r="workflow/scripts/functions/seuratDE.R",
  output:
    deseq2_pfvliver_list_rds="output/q1_pf_characterization/analyses/hc_pfvliver_manual_l3_deseq2_list.Rds",
    degs_pfvliver_xlsx="output/q1_pf_characterization/analyses/hc_pfvliver_manual_l3_degs.xlsx",
    deseq2_pfvcolon_list_rds="output/q1_pf_characterization/analyses/hc_pfvcolon_manual_l3_deseq2_list.Rds",
    degs_pfvcolon_xlsx="output/q1_pf_characterization/analyses/hc_pfvcolon_manual_l3_degs.xlsx",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/analyses/de_hc_pfvlivervcolon_manual_l3.log",
  benchmark:
    "output/q1_pf_characterization/analyses/de_hc_hc_pfvlivervcolon_manual_l3_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/de_hc_pfvlivervcolon.R "{input.pbmc_pf_liver_colon_macrophages_seuratobject_rds}" "{input.seuratDE_r}" "{output.deseq2_pfvliver_list_rds}" "{output.degs_pfvliver_xlsx}" "{output.deseq2_pfvcolon_list_rds}" "{output.degs_pfvcolon_xlsx}" &> "{log}"
    """

# rule de_hc_pf_m1vm2:
#   input:
#     hc_pf_mnp_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_mnp_SeuratObject.Rds",
#     seuratDE_r="workflow/scripts/functions/seuratDE.R",
#   output:
#     deseq2_list_rds="output/q1_pf_characterization/analyses/hc_pf_m1vm2_deseq2_list.Rds",
#   threads: 
#     1
#   conda:
#     "../envs/r-deseq2.yaml",
#   log:
#     "output/q1_pf_characterization/analyses/de_hc_pf_m1vm2.log",
#   benchmark:
#     "output/q1_pf_characterization/analyses/de_hc_pf_m1vm2_benchmark.txt",
#   resources:
#     mem_mb=60000,
#   shell:
#     """
#     Rscript --vanilla workflow/scripts/q1_pf_characterization/de_hc_pf_m1vm2.R "{input.hc_pf_mnp_seuratobject_rds}" "{input.seuratDE_r}" "{output.deseq2_list_rds}" &> "{log}"
#     """

###########
# Figures #
###########

rule fig_umap_hc_pbmc_pf_coltissue:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    umap_hc_pbmc_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_pbmc_pf_coltissue.pdf",
    umap_hc_pbmc_pf_coltissue_splittissue_pdf="output/q1_pf_characterization/figures/umap_hc_pbmc_pf_coltissue_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/umap_hc_pbmc_pf_coltissue.log",
  benchmark:
    "output/q1_pf_characterization/figures/umap_hc_pbmc_pf_coltissue_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_umap_hc_pbmc_pf_coltissue.R "{input.hc_pbmc_pf_seuratobject_rds}" "{output.umap_hc_pbmc_pf_coltissue_pdf}" "{output.umap_hc_pbmc_pf_coltissue_splittissue_pdf}" &> "{log}"
    """

rule fig_umap_hc_pbmc_pf_colcelltype_splittissue:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    umap_hc_pbmc_pf_colcelltype_splittissue_pdf="output/q1_pf_characterization/figures/umap_hc_pbmc_pf_col{level}_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/umap_hc_pbmc_pf_col{level}_splittissue.log",
  benchmark:
    "output/q1_pf_characterization/figures/umap_hc_pbmc_pf_col{level}_splittissue_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_umap_hc_pbmc_pf_colcelltype_splittissue.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.umap_hc_pbmc_pf_colcelltype_splittissue_pdf}" &> "{log}"
    """

rule fig_umap_hc_pf_colcelltype:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_gexumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_gex_pf_col{level}.pdf",
    hc_citeumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_cite_pf_col{level}.pdf",
    hc_wnnumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_wnn_pf_col{level}.pdf",
    hc_totalumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_total_pf_col{level}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/umap_hc_pf_col{level}.log",
  benchmark:
    "output/q1_pf_characterization/figures/umap_hc_pf_col{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_umap_hc_pf_colcelltype.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_gexumap_pf_coltissue_pdf}" "{output.hc_citeumap_pf_coltissue_pdf}" "{output.hc_wnnumap_pf_coltissue_pdf}" "{output.hc_totalumap_pf_coltissue_pdf}" &> "{log}"
    """
    
# rule fig_umap_hc_pf_splitmodality:
#   input:
#     hc_pf_subset_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_{subset}_SeuratObject.Rds",
#     celltype_markers_xlsx=config['celltype_markers'],
#   output:
#     hc_gexumap_pf_splitmodality_pdf="output/q1_pf_characterization/figures/umap_hc_gex_pf_col{level}.pdf",
#     hc_citeumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_cite_pf_col{level}.pdf",
#     hc_wnnumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_wnn_pf_col{level}.pdf",
#     hc_totalumap_pf_coltissue_pdf="output/q1_pf_characterization/figures/umap_hc_total_pf_col{level}.pdf",
#   threads: 
#     1
#   conda:
#     "../envs/r.yaml",
#   log:
#     "output/q1_pf_characterization/figures/umap_hc_pf_col{level}.log",
#   benchmark:
#     "output/q1_pf_characterization/figures/umap_hc_pf_col{level}_benchmark.txt",
#   resources:
#     mem_mb=16000,
#   params:
#     level="{level}",
#   shell:
#     """
#     Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_umap_hc_pf_colcelltype.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_gexumap_pf_coltissue_pdf}" "{output.hc_citeumap_pf_coltissue_pdf}" "{output.hc_wnnumap_pf_coltissue_pdf}" "{output.hc_totalumap_pf_coltissue_pdf}" &> "{log}"
#     """

rule fig_dotplot_markerexpression_hc_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    dotplot_markerexpression_hc_pbmc_pf_pdf="output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pbmc_pf_{level}.pdf",
    dotplot_markerexpression_hc_pbmc_pf_splittissue_pdf="output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pbmc_pf_{level}_splittissue.pdf",
    dotplot_markerexpression_hc_pf_gex_pdf="output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pf_{level}_gex.pdf",
    dotplot_markerexpression_hc_pf_pex_pdf="output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pf_{level}_pex.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pbmc_pf_{level}.log",
  benchmark:
    "output/q1_pf_characterization/figures/dotplot_markerexpression_hc_pbmc_pf_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_dotplot_markerexpression_hc_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.dotplot_markerexpression_hc_pbmc_pf_pdf}" "{output.dotplot_markerexpression_hc_pbmc_pf_splittissue_pdf}" "{output.dotplot_markerexpression_hc_pf_gex_pdf}" "{output.dotplot_markerexpression_hc_pf_pex_pdf}" &> "{log}"
    """

rule fig_heatmap_markerexpression_hc_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    heatmap_markerexpression_hc_pbmc_pf_pdf="output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pbmc_pf_{level}.pdf",
    heatmap_markerexpression_hc_pbmc_pf_splittissue_pdf="output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pbmc_pf_{level}_splittissue.pdf",
    heatmap_markerexpression_hc_pf_gex_pdf="output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pf_{level}_gex.pdf",
    heatmap_markerexpression_hc_pf_pex_pdf="output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pf_{level}_pex.pdf",
    heatmap_markerexpression_hc_pf_gex_pex_pdf="output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pf_{level}_gex_pex.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pbmc_pf_{level}.log",
  benchmark:
    "output/q1_pf_characterization/figures/heatmap_markerexpression_hc_pbmc_pf_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_heatmap_markerexpression_hc_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.heatmap_markerexpression_hc_pbmc_pf_pdf}" "{output.heatmap_markerexpression_hc_pbmc_pf_splittissue_pdf}" "{output.heatmap_markerexpression_hc_pf_gex_pdf}" "{output.heatmap_markerexpression_hc_pf_pex_pdf}" "{output.heatmap_markerexpression_hc_pf_gex_pex_pdf}" &> "{log}"
    """

rule fig_violinplot_hc_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    violinplot_hc_pbmc_pf_splittissue_pdf="output/q1_pf_characterization/figures/violinplot_hc_pbmc_pf_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/violinplot_hc_pbmc_pf.log",
  benchmark:
    "output/q1_pf_characterization/figures/violinplot_hc_pbmc_pf_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_violinplot_hc_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.violinplot_hc_pbmc_pf_splittissue_pdf}" &> "{log}"
    """

rule fig_hc_pfvpbmc_de:
  input:
    deseq2_list_rds="output/q1_pf_characterization/analyses/hc_pfvpbmc_{level}_deseq2_list.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    arrowplot_top5sig_pdf="output/q1_pf_characterization/figures/arrowplot_hc_pfvpbmc_{level}_top5sig.pdf",
    arrowplot_resident_pdf="output/q1_pf_characterization/figures/arrowplot_hc_pfvpbmc_{level}_resident.pdf",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/q1_pf_characterization/figures/arrowplot_hc_pfvpbmc_{level}_de.log",
  benchmark:
    "output/q1_pf_characterization/figures/arrowplot_hc_pfvpbmc_{level}_de_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_hc_pfvpbmc_de.R "{input.deseq2_list_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.arrowplot_top5sig_pdf}" "{output.arrowplot_resident_pdf}" &> "{log}"
    """

rule fig_boxplot_hc_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    boxplot_hc_pbmc_abundance_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_abundance_{comparison}.pdf",
    boxplot_hc_pbmc_abundance_patanno_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_abundance_patanno_{comparison}.pdf",
    boxplot_hc_pf_abundance_pdf="output/q1_pf_characterization/figures/boxplot_hc_pf_abundance_{comparison}.pdf",
    boxplot_hc_pf_abundance_patanno_pdf="output/q1_pf_characterization/figures/boxplot_hc_pf_abundance_patanno_{comparison}.pdf",
    boxplot_hc_pbmc_pf_abundance_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_abundance_{comparison}.pdf",
    boxplot_hc_pbmc_pf_abundance_patanno_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_abundance_patanno_{comparison}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_boxplot_hc_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.boxplot_hc_pbmc_abundance_pdf}" "{output.boxplot_hc_pbmc_abundance_patanno_pdf}" "{output.boxplot_hc_pf_abundance_pdf}" "{output.boxplot_hc_pf_abundance_patanno_pdf}" "{output.boxplot_hc_pbmc_pf_abundance_pdf}" "{output.boxplot_hc_pbmc_pf_abundance_patanno_pdf}" &> "{log}"
    """

rule fig_boxplot_hc_pbmc_pf_liver_colon_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    hc_liver_seuratobject_rds="resources/liver/liver_seuratobject.rds",
    hc_colon_seuratobject_rds="resources/colon/colon_seuratobject.rds",
    marker_genes_xlsx=config['celltype_markers'],
  output:
    boxplot_hc_pbmc_pf_colon_liver_abundance_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_colon_liver_abundance.pdf",
    boxplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf="output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_liver_colon_abundance.log",
  benchmark:
    "output/q1_pf_characterization/figures/boxplot_hc_pbmc_pf_liver_colon_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_boxplot_hc_pbmc_pf_liver_colon_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.hc_liver_seuratobject_rds}" "{input.hc_colon_seuratobject_rds}" "{input.marker_genes_xlsx}" "{output.boxplot_hc_pbmc_pf_colon_liver_abundance_pdf}" "{output.boxplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf}" &> "{log}"
    """

rule fig_stackedbarplot_hc_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    stackedbarplot_hc_pbmc_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_abundance_{level}.pdf",
    stackedbarplot_hc_pf_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pf_abundance_{level}.pdf",
    stackedbarplot_hc_pbmc_pf_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_abundance_{level}.pdf",
    stackedbarplot_hc_pbmc_abundance_myeloid_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_abundance_{level}_myeloid.pdf",
    stackedbarplot_hc_pf_abundance_myeloid_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pf_abundance_{level}_myeloid.pdf",
    stackedbarplot_hc_pbmc_pf_abundance_myeloid_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_abundance_{level}_myeloid.pdf",
    stackedbarplot_hc_pbmc_abundance_t_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_abundance_{level}_t.pdf",
    stackedbarplot_hc_pf_abundance_t_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pf_abundance_{level}_t.pdf",
    stackedbarplot_hc_pbmc_pf_abundance_t_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_abundance_{level}_t.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_abundance_{level}.log",
  benchmark:
    "output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_stackedbarplot_hc_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.level}" "{output.stackedbarplot_hc_pbmc_abundance_pdf}" "{output.stackedbarplot_hc_pf_abundance_pdf}" "{output.stackedbarplot_hc_pbmc_pf_abundance_pdf}" "{output.stackedbarplot_hc_pbmc_abundance_myeloid_pdf}" "{output.stackedbarplot_hc_pf_abundance_myeloid_pdf}" "{output.stackedbarplot_hc_pbmc_pf_abundance_myeloid_pdf}" "{output.stackedbarplot_hc_pbmc_abundance_t_pdf}" "{output.stackedbarplot_hc_pf_abundance_t_pdf}" "{output.stackedbarplot_hc_pbmc_pf_abundance_t_pdf}" &> "{log}"
    """

rule fig_stackedbarplot_hc_pbmc_pf_liver_colon_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
    hc_liver_seuratobject_rds="resources/liver/liver_seuratobject.rds",
    hc_colon_seuratobject_rds="resources/colon/colon_seuratobject.rds",
    marker_genes_xlsx=config['celltype_markers'],
  output:
    stackedbarplot_hc_pbmc_pf_colon_liver_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_colon_liver_abundance.pdf",
    stackedbarplot_hc_pbmc_pf_colon_woB_liver_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_colon_woB_liver_abundance.pdf",
    stackedbarplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf="output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_liver_colon_abundance.log",
  benchmark:
    "output/q1_pf_characterization/figures/stackedbarplot_hc_pbmc_pf_liver_colon_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_stackedbarplot_hc_pbmc_pf_liver_colon_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.hc_liver_seuratobject_rds}" "{input.hc_colon_seuratobject_rds}" "{input.marker_genes_xlsx}" "{output.stackedbarplot_hc_pbmc_pf_colon_liver_abundance_pdf}" "{output.stackedbarplot_hc_pbmc_pf_colon_woB_liver_abundance_pdf}" "{output.stackedbarplot_hc_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf}" &> "{log}"
    """

rule fig_heatmap_hc_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    heatmapplot_hc_pbmc_pf_abundance_scaled_pdf="output/q1_pf_characterization/figures/heatmap_hc_pbmc_pf_abundance_{comparison}_scaled.pdf",
    heatmapplot_hc_pbmc_pf_abundance_unscaled_pdf="output/q1_pf_characterization/figures/heatmap_hc_pbmc_pf_abundance_{comparison}_unscaled.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/q1_pf_characterization/figures/heatmap_hc_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/q1_pf_characterization/figures/heatmap_hc_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_heatmap_hc_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.heatmapplot_hc_pbmc_pf_abundance_scaled_pdf}" "{output.heatmapplot_hc_pbmc_pf_abundance_unscaled_pdf}" &> "{log}"
    """

# rule fig_hc_mnp_markerexpression_pf:
#   input:
#     hc_pf_mnp_seuratobject_rds="output/q1_pf_characterization/subsets/hc_pf_MNP_SeuratObject.Rds",
#   output:
#     violinplot_hc_mnp_markerexpression_pf_pdf="output/q1_pf_characterization/figures/violinplot_hc_mnp_markerexpression_pf.pdf",
#     boxplot_hc_mnp_markerexpression_pf_pdf="output/q1_pf_characterization/figures/boxplot_hc_mnp_markerexpression_pf.pdf",
#     dotplot_hc_mnp_markerexpression_pf_pdf="output/q1_pf_characterization/figures/dotplot_hc_mnp_markerexpression_pf.pdf",
#     umap_hc_mnp_markerexpression_pf_pdf="output/q1_pf_characterization/figures/umap_hc_mnp_markerexpression_pf.pdf",
#   threads: 
#     1
#   conda:
#     "../envs/r.yaml",
#   log:
#     "output/q1_pf_characterization/figures/hc_mnp_markerexpression_pf.log",
#   benchmark:
#     "output/q1_pf_characterization/figures/hc_mnp_markerexpression_pf_benchmark.txt",
#   resources:
#     mem_mb=16000,
#   shell:
#     """
#     Rscript --vanilla workflow/scripts/q1_pf_characterization/fig_hc_mnp_markerexpression_pf.R "{input.hc_pf_mnp_seuratobject_rds}" "{output.violinplot_hc_mnp_markerexpression_pf_pdf}" "{output.boxplot_hc_mnp_markerexpression_pf_pdf}" "{output.dotplot_hc_mnp_markerexpression_pf_pdf}" "{output.umap_hc_mnp_markerexpression_pf_pdf}" &> "{log}"
#     """
