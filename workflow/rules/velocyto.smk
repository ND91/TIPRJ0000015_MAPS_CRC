rule velocyto:
  input:
    cellrangerdir="output/cellranger/{run}",
    genesgtf="resources/reference_genome/genes/genes.gtf",
    repeatmaskgtf=config['repeatmask_file'],
  output:
    loomfile="output/cellranger/{run}/velocyto/{run}.loom",
  conda:
    "../envs/velocyto.yaml",
  log:
    "output/cellranger/{run}_velocyto.log",
  threads: 
    8
  benchmark:
    "output/cellranger/{run}_velocyto_benchmark.txt",
  resources:
    mem_mb=64000,
  shell:
    """
    velocyto run10x -m "{repeatmaskgtf}" "{cellrangerdir}" "{genesgtf}"
    """