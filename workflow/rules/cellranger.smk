rule cellranger_count:
  input:
    fastqdir="resources/fastq/{run}",
    reference_genome="resources/reference_genome",
  output:
    cellrangerdir=touch(directory("output/cellranger/{run}")),
    filtered_matrix=touch(directory("output/cellranger/{run}/outs/filtered_feature_bc_matrix")),
    raw_matrix=touch(directory("output/cellranger/{run}/outs/raw_feature_bc_matrix")),
  log:
    "output/cellranger/{run}_cellranger.log",
  threads: 
    8
  benchmark:
    "output/cellranger/{run}_count.benchmark.txt"
  params:
    runid="{run}",
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