rule velocyto:
  input:
    #bc_tsv="output/velocyto/{run}/filtered_barcodes.tsv",
    cellrangerdir="output/cellranger/{run}",
    # bam_cellpossorted="output/velocyto/{run}/cellsorted_possorted_genome_bam.bam",
    # #bai_cellpossorted="output/velocyto/{run}/cellsorted_possorted_genome_bam.bam.bai",
    genesgtf="resources/reference_genome/genes/genes.gtf",
    repeatmaskgtf=config['repeatmask_file'],
  output:
    loomfile="output/cellranger/{run}/velocyto/{run}.loom",
  conda:
    "../envs/velocyto.yaml",
  log:
    "output/cellranger/{run}/velocyto.log",
  threads: 
    8
  benchmark:
    "output/cellranger/{run}/velocyto_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    velocyto run10x --samtools-threads {threads} --samtools-memory 8000 -m "{input.repeatmaskgtf}" "{input.cellrangerdir}" "{input.genesgtf}" &> "{log}"
    """

rule velocyto_merge_loom:
  input:
    looms=expand("output/cellranger/{run}/velocyto/{run}.loom", run=runs),
  output:
    velocyto_merged_loom="output/merged/velocyto_merged.loom",
  conda:
    "../envs/velocyto.yaml",
  log:
    "output/merged/merge_loom.log",
  threads: 
    8
  benchmark:
    "output/merged/merge_loom_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    python3 workflow/scripts/preparation/merge_loom.py --velocyto_looms {input.looms:q} --combined_loom {output.velocyto_merged_loom:q} &> {log:q}
    """

rule loom_curate:
  input:
    velocyto_merged_loom="output/merged/velocyto_merged.loom",
    metadata_csv="output/curated/SeuratObject_metadata.csv",
  output:
    velocyto_curated_loom="output/curated/velocyto_curated.loom",
  conda:
    "../envs/velocyto.yaml",
  log:
    "output/curated/loom_curate.log",
  threads: 
    8
  benchmark:
    "output/merged/loom_curate_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    python3 workflow/scripts/preparation/loom_curate.py --merged_loom "{input.velocyto_merged_loom}" --cellmetadata_csv "{input.metadata_csv}" --annotated_loom "{output.velocyto_curated_loom}" &> "{log}"
    """

# rule loom_to_sce:
#   input:
#     velocyto_merged_loom="output/merged/velocyto_merged.loom",
#   output:
#     velocyto_sce_rds="output/merged/velocyto_merged_sce.Rds",
#   conda:
#     "../envs/velocyto.yaml",
#   log:
#     "output/merged/loom_to_sce.log",
#   threads: 
#     8
#   benchmark:
#     "output/merged/loom_to_sce_benchmark.txt",
#   resources:
#     mem_mb=64000,
#   shell:
#     """
#     Rscript workflow/scripts/preparation/loom_to_sce.R "{input.velocyto_merged_loom:q}" "{output.velocyto_sce_rds:q}" &> "{log}"
#     """
