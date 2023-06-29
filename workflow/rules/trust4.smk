rule trust4:
  input:
    cellrangerdir="output/cellranger/{run}",
    hg38_bcrtcr_fa=config['bcrtcr'],
  output:
    trust4dir=touch(directory("output/trust4/{run}")),
    airr_tsv = "output/trust4/{run}/{run}_barcode_airr.tsv",
    # annot_fa="output/trust4/{run}/{run}_annot.fa",
    # assembled_reads_fa="output/trust4/{run}/{run}_assembled_reads.fa",
    # barcode_report_tsv="output/trust4/{run}/{run}_barcode_report.tsv",
    # cd3_out="output/trust4/{run}/{run}_cdr3.out",
    # final_out="output/trust4/{run}/{run}_final.out",
    # raw_out="output/trust4/{run}/{run}_raw.out",
    # report_tsv="output/trust4/{run}/{run}_report.tsv",
    # toassemble_fq="output/trust4/{run}/{run}_toassemble.fq",
    # toassemble_bc_fa="output/trust4/{run}/{run}_toassemble_bc.fa",
  conda:
    "../envs/trust4.yaml",
  log:
    "output/trust4/{run}/trust4.log",
  threads: 
    8
  benchmark:
    "output/trust4/{run}/trust4_benchmark.txt",
  resources:
    mem_mb=10000,
  params:
    run="{run}",
  shell:
    """
    run-trust4 -b "{input.cellrangerdir}/outs/possorted_genome_bam.bam" --barcode CB -t {threads} -f "{input.hg38_bcrtcr_fa}" -o "{params.run}" --od "{output.trust4dir}" &> "{log}"
    """

rule trust4_annotate:
  input:
    airr_tsv = "output/trust4/{run}/{run}_barcode_airr.tsv",
  output:
    airr_annotated_csv = "output/trust4/{run}/{run}_barcode_annotated_airr.csv",
  conda:
    "../envs/r.yaml",
  log:
    "output/trust4/{run}/trust4_annotate.log",
  threads: 
    1
  benchmark:
    "output/trust4/{run}/trust4_annotate_benchmark.txt",
  resources:
    mem_mb=10000,
  params:
    runid="{run}"
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/trust4_annotate.R "{input.airr_tsv}" "{params.runid}" "{output.airr_annotated_csv}" &> "{log}"
    """

rule trust4_merge:
  input:
    airr_annotated_csvs=expand("output/trust4/{run}/{run}_barcode_annotated_airr.csv", run=runs), 
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    airr_merged_csv="output/trust4/trust4_barcode_annotated_airr.csv",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/trust4/trust4_merge.log",
  benchmark:
    "output/trust4/trust4_merge_benchmark.txt",
  resources:
    mem_mb=30000,
  shell:
    """
    Rscript --vanilla workflow/scripts/preparation/trust4_merge.R "{output.airr_merged_csv}" "{input.live_singlet_nonproliferating_seuratobject_rds}" {input.airr_annotated_csvs} &> "{log}"
    """
