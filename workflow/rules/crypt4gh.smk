rule crypt4gh:
  input:
    fastqdir=lambda w: (fastq_df[fastq_df.file == w.fastq_fn].old_path).tolist()
  output:
    fastqc4ghdir=touch(directory("output/cellranger/{run_nonfbc}")),
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
