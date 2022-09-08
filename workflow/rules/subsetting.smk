# Tissues: PBMC, PF, TX

rule pbmc_pf_tx_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/pbmc_pf_tx/pbmc_pf_tx_subsetting.log",
  benchmark:
    "output/pbmc_pf_tx/pbmc_pf_tx_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_seuratobject_rds} &> {log}
    """

# Tissues: PBMC, PF, TX. Cells: alive, non-debris, non-doublets, non-proliferating

rule pbmc_pf_tx_cleaned_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_cleaned_seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/pbmc_pf_tx_cleaned_subsetting.log",
  benchmark:
    "output/subsets/pbmc_pf_tx_cleaned_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_cleaned_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_cleaned_seuratobject_rds} &> {log}
    """

# Tissues: PBMC, PF, TX. Cells: CD45+, alive, non-debris, non-doublets, non-proliferating

rule pbmc_pf_tx_cleaned_cd45p_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_cleaned_cd45p_seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/pbmc_pf_tx_cleaned_cd45p_subsetting.log",
  benchmark:
    "output/subsets/pbmc_pf_tx_cleaned_cd45p_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_cleaned_cd45p_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_cleaned_cd45p_seuratobject_rds} &> {log}
    """

# Tissues: PBMC, PF. Patients: non-metastatic. Cells: CD45+, alive, non-debris, non-doublets, non-proliferating

rule pbmc_pf_tx_cleaned_cd45p_pmn_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_cleaned_cd45p_pmn_seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/pbmc_pf_tx_cleaned_cd45p_pmn_subsetting.log",
  benchmark:
    "output/subsets/pbmc_pf_tx_cleaned_cd45p_pmn_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_cleaned_cd45p_pmn_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_cleaned_cd45p_pmn_seuratobject_rds} &> {log}
    """

# Tissues: PBMC, PF, TX. Patients: Common donors 

rule pbmc_pf_tx_paired_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_paired_seuratobject_rds="output/subsets/pbmc_pf_tx_paired_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/pbmc_pf_tx_paired_subsetting.log",
  benchmark:
    "output/subsets/pbmc_pf_tx_paired_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_paired_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_paired_seuratobject_rds} &> {log}
    """

# Tissues: PBMC, PF, TX. Patients: Common donors. Cells: CD45+, alive, non-debris, non-doublets, non-proliferating

rule pbmc_pf_tx_paired_cleaned_cd45p_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    pbmc_pf_tx_paired_cleaned_cd45p_seuratobject_rds="output/subsets/pbmc_pf_tx_paired_cleaned_cd45p_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/pbmc_pf_tx_paired_cleaned_cd45p_subsetting.log",
  benchmark:
    "output/subsets/pbmc_pf_tx_paired_cleaned_cd45p_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/pbmc_pf_tx_paired_cleaned_cd45p_subsetting.R {input.seurat_curated_rds} {output.pbmc_pf_tx_paired_cleaned_cd45p_seuratobject_rds} &> {log}
    """

# Tissues: TX 

rule tx_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    tx_seuratobject_rds="output/subsets/tx_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/tx_subsetting.log",
  benchmark:
    "output/subsets/tx_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/tx_subsetting.R {input.seurat_curated_rds} {output.tx_seuratobject_rds} &> {log}
    """

# Tissues: TX. Cells: CD45+, alive, non-debris, non-doublets, non-proliferating

rule tx_cleaned_cd45p_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    tx_cleaned_cd45p_seuratobject_rds="output/subsets/tx_cleaned_cd45p_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/tx_cleaned_cd45p_subsetting.log",
  benchmark:
    "output/subsets/tx_cleaned_cd45p_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/tx_cleaned_cd45p_subsetting.R {input.seurat_curated_rds} {output.tx_cleaned_cd45p_seuratobject_rds} &> {log}
    """

# Cells: CD4T

rule cd4t_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    cd4t_seuratobject_rds="output/subsets/cd4t_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/cd4t_subsetting.log",
  benchmark:
    "output/subsets/cd4t_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/subsetting/cd4t_subsetting.R {input.seurat_curated_rds} {output.tx_seuratobject_rds} &> {log}
    """