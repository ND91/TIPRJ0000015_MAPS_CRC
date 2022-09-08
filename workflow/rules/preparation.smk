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

rule celltype_curation:
  input:
    seurat_merged_rds="output/merged/merged_SeuratObject.Rds",
    curated_csv=config["curated_celltypes"],
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