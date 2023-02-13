rule cellranger_count_nonfbc:
  input:
    fastqdir="resources/fastq/{run_nonfbc}",
    reference_genome="resources/reference_genome",
  output:
    cellrangerdir=touch(directory("output/cellranger/{run_nonfbc}")),
    filtered_matrix=touch(directory("output/cellranger/{run_nonfbc}/outs/filtered_feature_bc_matrix")),
    raw_matrix=touch(directory("output/cellranger/{run_nonfbc}/outs/raw_feature_bc_matrix")),
    bc_tsv_gz="output/cellranger/{run_nonfbc}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    bam_possorted="output/cellranger/{run_nonfbc}/outs/possorted_genome_bam.bam",
  log:
    "output/cellranger/{run_nonfbc}_cellranger.log",
  threads: 
    8
  benchmark:
    "output/cellranger/{run_nonfbc}_count.benchmark.txt"
  params:
    runid="{run_nonfbc}",
    cellranger_basepath=config["cellranger_basepath"],
  resources:
    mem_mb=64000,
  shell:
    """
    rm -rf {output.cellrangerdir}
    mem_gb=$(echo {resources.mem_mb}/1024|bc)
    mkdir -p "output/cellranger"
    cd "output/cellranger"
    #find . -type d -empty -delete
    {params.cellranger_basepath}/cellranger count \
    --id="{params.runid}" \
    --fastqs=$OLDPWD/{input.fastqdir} \
    --transcriptome=$OLDPWD/{input.reference_genome} \
    --nosecondary \
    --disable-ui \
    --localcores={threads} \
    --localmem=$mem_gb &> $OLDPWD/{log}
    """

rule cellranger_count_fbc:
  input:
    #fastqdir=lambda wildcards: 'resources/fastq/' + lookup_fbc.loc[lookup_fbc['Run'] == wildcards.run_fbc, "Run_well"].astype(str),
    #fastqfbcdir=lambda wildcards: 'resources/fastq/' + lookup_fbc.loc[lookup_fbc['Run'] == wildcards.run_fbc, "Run_well"].astype(str) + '-ADT',
    fastqdir="resources/fastq/{run_fbc}",
    fastqfbcdir="resources/fastq/{run_fbc}-ADT",
    libraries="config/samples/fbcs/{run_fbc}_libraries.csv",
    features="config/samples/fbcs/TotalSeq_B_Human_Universal_Cocktail_V1_399904_Antibody_reference_UMI_counting_CellRanger.csv",
    reference_genome="resources/reference_genome",
  output:
    cellrangerdir=touch(directory("output/cellranger/{run_fbc}")),
    filtered_matrix=touch(directory("output/cellranger/{run_fbc}/outs/filtered_feature_bc_matrix")),
    raw_matrix=touch(directory("output/cellranger/{run_fbc}/outs/raw_feature_bc_matrix")),
    bc_tsv_gz="output/cellranger/{run_fbc}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    bam_possorted="output/cellranger/{run_fbc}/outs/possorted_genome_bam.bam",
  log:
    "output/cellranger/{run_fbc}_cellranger.log",
  threads: 
    8
  benchmark:
    "output/cellranger/{run_fbc}_count.benchmark.txt"
  params:
    runid="{run_fbc}",
    cellranger_basepath=config["cellranger_basepath"],
  resources:
    mem_mb=64000,
  shell:
    """
    rm -rf {output.cellrangerdir}
    mem_gb=$(echo {resources.mem_mb}/1024|bc)
    mkdir -p {output.cellrangerdir}
    cd "output/cellranger"
    #find . -type d -empty -delete
    {params.cellranger_basepath}/cellranger count \
      --id="{params.runid}" \
      --libraries=$OLDPWD/{input.libraries} \
      --feature-ref=$OLDPWD/{input.features} \
      --transcriptome=$OLDPWD/{input.reference_genome} \
      --nosecondary \
      --disable-ui \
      --localcores={threads} \
      --localmem=$mem_gb &> $OLDPWD/{log}
    """
