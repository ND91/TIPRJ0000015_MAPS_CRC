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

rule download_gca_colon:
  output:
    colon_gca_h5ad="resources/reference_data/colon_gca.h5ad",
  threads:
    1
  conda:
    "../envs/wget.yaml",
  message:
    "--- Downloading Gut Cell Atlas Colon Immune Atlas ---"
  log:
    "resources/reference_data/download_colon_gca_h5ad.log",
  benchmark:
    "resources/reference_data/download_colon_gcat_benchmark.txt",
  shell:
    """
    wget --continue -O {output.colon_gca_h5ad} 'https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Colon_cell_atlas.h5ad'
    """

rule download_hca_pan_immune_system:
  output:
    pan_immune_hca_h5ad="resources/reference_data/pan_immune_gca.h5ad",
  threads:
    1
  conda:
    "../envs/wget.yaml",
  message:
    "--- Downloading Human Cell Atlas Pan Immune System Atlas ---"
  log:
    "resources/reference_data/download_hca_pan_immune_system_h5ad.log",
  benchmark:
    "resources/reference_data/download_hca_pan_immune_system_h5ad_benchmark.txt",
  shell:
    """
    wget --continue -O {output.pan_immune_hca_h5ad} 'https://storage.googleapis.com/datarepo-be8a010d-bucket/afb53e46-8262-4a81-a1c7-8f93ee63eede/7a9fccf1-764d-42e0-a628-4daae8fbb229/PAN.A01.v01.raw_count.20210429.PFI.embedding.h5ad?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230213%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230213T122753Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-17b90eb5&X-Goog-Signature=9569308c31cca32843b1bd208bc0714271599b22f5d4ad8826833a1ff70afd309522904053f653f90d0006d1d8bd9eb3cfa29cfa7dc39a89bdb250b296c9acf34b74ddb197662cf4dde9512f5f79135165644a4b9fb10d11d0a16f85b70c777efe6881858e65050bccc33e3105bd5685cd526bb06958de2fc991b3969bdd082421a497f4414c5e74a0022c84a2b8294595e1a2e19c68a871eb8cecd0ce7839f265aecbc7382486d74166d30a41024d2c529004bf37a6cd177c0783bde1a97c5370b2fa0b7be239259293f91f753067c97c270cd7fdc8dd7efa30d7f40d2cbd3934bf1a71530be3d97882d94bf6801a686bb40e76b021bf99d9f7e56e926733b7'
    """

rule gse134355_preparation:
  input:
    gse134355_seuratobject_rds=config['gse134355_seuratobject_rds'],
    gse134355_celltranslations_xlsx=config['gse134355_celltranslations_xlsx'],
  output:
    gse134355_reannotated_seurat_rds="resources/GSE134355/gse134355_reannotated_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "output/resources/GSE134355/gse134355_preparation.log",
  benchmark:
    "output/resources/GSE134355/gse134355_preparation_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE134355: Preparing ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse134355_preparation.R "{input.gse134355_seuratobject_rds}" "{input.gse134355_celltranslations_xlsx}" "{output.gse134355_reannotated_seurat_rds}" &> "{log}"
    """

rule gse134355_subsetting_adult_immune:
  input:
    gse134355_reannotated_seurat_rds="resources/GSE134355/gse134355_reannotated_SeuratObject.Rds",
  output:
    gse134355_adult_immune_seurat_rds="resources/GSE134355/gse134355_adult_immune_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "resources/GSE134355/gse134355_subsetting_adult_immune.log",
  benchmark:
    "resources/GSE134355/gse134355_subsetting_adult_immune_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE134355: Subsetting adult immune ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse134355_subset_adult_immune.R "{input.gse134355_reannotated_seurat_rds}" "{output.gse134355_adult_immune_seurat_rds}" &> "{log}"
    """

rule gse134355_adult_immune_celltype_curation:
  input:
    gse134355_adult_immune_seurat_rds="resources/GSE134355/gse134355_adult_immune_SeuratObject.Rds",
    gse134355_curated_csv=config["gse134355_curated_adult_immune_celltypes"],
  output:
    gse134355_adult_immune_reannotated_seurat_rds="resources/GSE134355/gse134355_adult_immune_reannotated_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/GSE134355/gse134355_adult_immune_celltype_curation.log",
  benchmark:
    "output/GSE134355/gse134355_adult_immune_celltype_curation_benchmark.txt",
  resources:
    mem_mb=40000,
  shell:
    """
    Rscript --vanilla workflow/scripts/resources/gse134355_adult_immune_celltype_curation.R "{input.gse134355_adult_immune_seurat_rds}" "{input.gse134355_curated_csv}" "{output.gse134355_adult_immune_reannotated_seurat_rds}" &> "{log}"
    """

rule gse178318_preparation:
  input:
    gse178318_mtx=config['gse178318_mtx'],
    gse178318_genes=config['gse178318_genes'],
    gse178318_barcodes=config['gse178318_barcodes'],
    gse178318_curated_celltypes=config['gse178318_curated_celltypes'],
  output:
    gse178318_annotated_seurat_rds="resources/GSE178318/gse178318_annotated_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "output/resources/GSE178318/gse178318_preparation.log",
  benchmark:
    "output/resources/GSE178318/gse178318_preparation_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE178318: Preparing ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse178318_preparation.R "{input.gse178318_mtx}" "{input.gse178318_genes}" "{input.gse178318_barcodes}" "{input.gse178318_curated_celltypes}" "{output.gse178318_annotated_seurat_rds}" &> "{log}"
    """

rule gse178318_subsetting_mnp:
  input:
    gse178318_annotated_seurat_rds="resources/GSE178318/gse178318_annotated_SeuratObject.Rds",
  output:
    gse178318_mnp_seurat_rds="resources/GSE178318/gse178318_mnp_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "output/resources/GSE178318/gse178318_subsetting_mnp.log",
  benchmark:
    "output/resources/GSE178318/gse178318_subsetting_mnp_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE178318: Subsetting MNP ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse178318_subset_mnp.R "{input.gse178318_annotated_seurat_rds}" "{output.gse178318_mnp_seurat_rds}" &> "{log}"
    """

rule gse178318_subsetting_macrophages:
  input:
    gse178318_annotated_seurat_rds="resources/GSE178318/gse178318_annotated_SeuratObject.Rds",
  output:
    gse178318_macrophages_seurat_rds="resources/GSE178318/gse178318_macrophages_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "output/resources/GSE178318/gse178318_subsetting_macrophages.log",
  benchmark:
    "output/resources/GSE178318/gse178318_subsetting_macrophages_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE178318: Subsetting macrophages ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse178318_subset_macrophages.R "{input.gse178318_annotated_seurat_rds}" "{output.gse178318_macrophages_seurat_rds}" &> "{log}"
    """

rule gse201333_subsetting_lung:
  input:
    gse201333_seurat_rds=config['gse201333_seuratobject_rds'],
  output:
    gse201333_lung_seurat_rds="resources/GSE201333/gse201333_lung_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "resources/GSE201333/gse201333_subsetting_lung.log",
  benchmark:
    "resources/GSE201333/gse201333_subsetting_lung_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE201333: Subsetting lung ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse201333_subset_lung.R "{input.gse201333_seurat_rds}" "{output.gse201333_lung_seurat_rds}" &> "{log}"
    """

rule gse201333_subsetting_lymphnode:
  input:
    gse201333_seurat_rds=config['gse201333_seuratobject_rds'],
  output:
    gse201333_lymphnode_seurat_rds="resources/GSE201333/gse201333_lymphnode_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "resources/GSE201333/gse201333_subsetting_lymphnode.log",
  benchmark:
    "resources/GSE201333/gse201333_subsetting_lymphnode_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE201333: Subsetting lymph node ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse201333_subset_lymphnode.R "{input.gse201333_seurat_rds}" "{output.gse201333_lymphnode_seurat_rds}" &> "{log}"
    """

rule gse201333_preparation:
  input:
    gse201333_seuratobject_rds=config['gse201333_seuratobject_rds'],
    gse201333_celltranslations_xlsx=config['gse201333_celltranslations_xlsx'],
  output:
    gse201333_reannotated_seurat_rds="resources/GSE201333/gse201333_reannotated_SeuratObject.Rds",
  conda:
    "../envs/r.yaml",
  log:
    "output/resources/GSE201333/gse201333_preparation.log",
  benchmark:
    "output/resources/GSE201333/gse201333_preparation_benchmark.txt",
  resources:
    mem_mb=60000,
  threads: 
    1
  message:
    "--- GSE201333: Preparing ---",
  shell:
    """
    Rscript workflow/scripts/resources/gse201333_preparation.R "{input.gse201333_seuratobject_rds}" "{input.gse201333_celltranslations_xlsx}" "{output.gse201333_reannotated_seurat_rds}" &> "{log}"
    """
