import os
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

rule copy_fastqs:
  input:
    fastqdir=lambda w: (fastq_df[fastq_df.file == w.fastq_fn].old_path).tolist()
  output:
    directory("resources/fastq/{fastq_fn}"),
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying {input.fastqdir} to scratch ---"
  shell:
    """
    if [ ! -d "{input.fastqdir}" ]; then
      echo "{input.fastqdir} does not exist!"
      exit 1
    fi
    
    rsync -ar {input.fastqdir}/ {output}
    """

rule copy_cellranger_reference:
  input:
    config['reference_dir']
  output:
    directory("resources/reference_genome")
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying reference genome to scratch ---"
  shell:
    """
    rsync -ar {input}/ {output}
    """

rule download_pbmc_reference_h5seurat:
  output:
    pbmc_reference_h5seurat="resources/reference_data/pbmc_multimodal.h5seurat",
  threads: 
    1
  conda:
    "../envs/wget.yaml",
  message:
    "--- Downloading PBMC reference dataset in h5seurat format from Hao et al. 2020 ---"
  log:
    "resources/reference_data/download_pbmc_reference_h5seurat.log",
  benchmark:
    "resources/reference_data/download_pbmc_reference_h5seurat_benchmark.txt",
  shell:
    """
    wget --continue -O {output.pbmc_reference_h5seurat} 'https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat'
    """