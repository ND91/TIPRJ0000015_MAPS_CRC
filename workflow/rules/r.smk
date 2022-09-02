rule preparation_pbmc_reference:
  input:
    pbmc_reference_h5seurat="resources/reference_data/pbmc_multimodal.h5seurat",
  output:
    pbmc_reference_rds="resources/reference_data/pbmc_multimodal.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml"
  log:
    "resources/reference_data/pbmc_reference_preparation.log",
  benchmark:
    "resources/reference_data/pbmc_reference_preparation_benchmark.txt",
  resources:
    mem_mb=40000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/pbmc_multimodal_preparation.R {input.pbmc_reference_h5seurat} {output.pbmc_reference_rds} &> {log}
    """

rule normalization:
  input:
    count_path="output/cellranger/{run}/outs/filtered_feature_bc_matrix",
  output:
    seurat_rds_path="output/normalized/{run}_normalized_SeuratObject.Rds",
  threads:
    1
  conda:
    "../envs/r.yaml"
  log:
    "output/normalized/{run}.log",
  benchmark:
    "output/normalized/{run}_benchmark.txt"
  resources:
    mem_mb=16000,
  params:
    runid="{run}",
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/normalization.R {input.count_path} {output.seurat_rds_path} {params.runid} &> {log}
    """

rule celltype_annotation:
  input:
    seurat_rds_path="output/normalized/{run}_normalized_SeuratObject.Rds",
    pbmc_rds_path="resources/pf_atlas/pf_atlas_seuratObject.Rds",
  output:
    seurat_annotated_rds_path="output/cell_metadata/celltype_annotation/{run}_celltype_annotated_SeuratObject.Rds",
    seurat_annotated_csv_path="output/cell_metadata/celltype_annotation/{run}_celltype_annotated.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/cell_metadata/celltype_annotation/{run}_celltype_annotation.log",
  benchmark:
    "output/cell_metadata/celltype_annotation/{run}_celltype_annotation_benchmark.txt",
  resources:
    mem_mb=47000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/celltype_annotation.R {input.seurat_rds_path} {input.pbmc_rds_path} {output.seurat_annotated_rds_path} {output.seurat_annotated_csv_path} &> {log}
    """

rule sample_metadata_annotate:
  input:
    seurat_nonfbc_rds="output/cell_metadata/celltype_annotation/{run_nonfbc}_celltype_annotated_SeuratObject.Rds",
    sample_metadata=config["sample_metadata"],
  output:
    seurat_sample_metadata_annotated_rds="output/cell_metadata/sample_metadata_annotation/{run_nonfbc}_sample_metadata_annotated_SeuratObject.Rds",
    seurat_sample_metadata_annotated_csv="output/cell_metadata/sample_metadata_annotation/{run_nonfbc}_sample_metadata_annotated.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/cell_metadata/sample_metadata_annotation/{run_nonfbc}_sample_metadata_annotated.log",
  benchmark:
    "output/cell_metadata/sample_metadata_annotation/{run_nonfbc}_sample_metadata_annotated_benchmark.txt"
  resources:
    mem_mb=12000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/nonfbc_sample_metadata_annotate.R {input.seurat_nonfbc_rds} {input.sample_metadata} {output.seurat_sample_metadata_annotated_rds} {output.seurat_sample_metadata_annotated_csv} &> {log}
    """
 
rule merge:
  input:
    seurat_rds=expand("output/cell_metadata/sample_metadata_annotation/{run}_sample_metadata_annotated_SeuratObject.Rds", run=runs), 
  output:
    seurat_merged_rds="output/merged/merged_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/merged/merge.log",
  benchmark:
    "output/merged/merge_benchmark.txt",
  resources:
    mem_mb=150000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/merge.R {output.seurat_merged_rds} {input.seurat_rds} &> {log}
    """
  
# Curations were done in a manual fashion where the automatic annotations were checked against canonical markers.

rule celltype_curation:
  input:
    seurat_merged_rds="output/merged/merged_SeuratObject.Rds",
    curated_csv=config['curated_celltypes']
  output:
    seurat_curated="output/curated/curated_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/curated/celltype_curation.log",
  benchmark:
    "output/curated/celltype_curation_benchmark.txt",
  resources:
    mem_mb=40000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/celltype_curation.R {input.seurat_merged_rds} {input.curated_csv} {output.seurat_curated} &> {log}
    """

# Analysis 1: TX composition

rule tx_preparation:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    tx_seuratobject_rds="output/tx/tx_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/tx/tx_preparation.log",
  benchmark:
    "output/tx/tx_preparation_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/tx_preparation.R {input.seurat_curated_rds} {output.tx_seuratobject_rds} &> {log}
    """
    
rule fig_boxplot_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/tx/tx_SeuratObject.Rds",
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
 
# Analysis 2: PBMC vs PF vs TX from paired patients 

rule pbmc_pf_tx_paired_preparation:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_paired_seuratobject_rds="output/pbmc_pf_tx_paired/pbmc_pf_tx_paired_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/pbmc_pf_tx_paired/pbmc_pf_tx_paired_preparation.log",
  benchmark:
    "output/pbmc_pf_tx_paired/pbmc_pf_tx_paired_preparation_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/pbmc_pf_tx_paired_preparation.R {input.seurat_curated_rds} {output.pbmc_pf_tx_paired_seuratobject_rds} &> {log}
    """
    
rule fig_pbmc_pf_tx_tsne:
  input:
    pbmc_pf_tx_paired_seuratobject_rds="output/pbmc_pf_tx_paired/pbmc_pf_tx_paired_SeuratObject.Rds",
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

rule txvpf_paired_de:
  input:
    pbmc_pf_tx_paired_seuratobject_rds="output/pbmc_pf_tx_paired/pbmc_pf_tx_paired_SeuratObject.Rds",
  output:
    degs_csv="output/pbmc_pf_tx_paired/{celltype}/{celltype}_degs.csv",
    deseq2_rds="output/pbmc_pf_tx_paired/{celltype}/{celltype}_deseq2.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/pbmc_pf_tx_paired/{celltype}/{celltype}_txvpf_paired_de.log",
  benchmark:
    "output/pbmc_pf_tx_paired/{celltype}/{celltype}_txvpf_paired_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    celltype_name="{celltype}",
    celltype_column_name="manual_l4",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_txvpf_paired_de.R {input.pbmc_pf_tx_paired_seuratobject_rds} "{params.celltype_column_name}" "{params.celltype_name}" "{output.degs_csv}" "{output.deseq2_rds}" &> "{log}"
    """    

rule fig_boxplot_pf_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/tx/tx_SeuratObject.Rds",
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

# Analysis 3: PBMC vs PF

rule fig_boxplot_pbmc_pf_l3rl2:
  input:
    seurat_curated="output/curated/curated_SeuratObject.Rds",
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
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_pbmc_pf_l3rl2.R "{input.seurat_curated}" "{output.boxplot_pbmc_pf_l3rl2_df_csv}" "{output.boxplot_pbmc_pf_l3rl2_svg}" &> "{log}"
    """
