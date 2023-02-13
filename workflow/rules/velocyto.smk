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
