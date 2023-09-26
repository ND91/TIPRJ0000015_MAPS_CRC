# CD45+, Living, singlet, non-proliferating cells

rule live_singlet_nonproliferating_subsetting:
  input:
    seurat_curated_rds="output/curated/curated_SeuratObject.Rds",
  output:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/live_singlet_nonproliferating_subsetting.log",
  benchmark:
    "output/subsets/live_singlet_nonproliferating_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/cd45p_live_singlet_nonproliferating_subsetting.R "{input.seurat_curated_rds}" "{output.live_singlet_nonproliferating_seuratobject_rds}" &> "{log}"
    """

## HC: PBMC and PF

rule hc_pbmc_pf_subsetting:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/hc_pbmc_pf_subsetting.log",
  benchmark:
    "output/subsets/hc_pbmc_pf_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/hc_pbmc_pf_subsetting.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_pbmc_pf_seuratobject_rds}" &> "{log}"
    """

rule hc_pbmc_pf_b_subsetting:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_pbmc_pf_b_seuratobject_rds="output/subsets/hc_pbmc_pf_B_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/hc_pbmc_pf_b_subsetting.log",
  benchmark:
    "output/subsets/hc_pbmc_pf_b_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/hc_pbmc_pf_b_subsetting.R "{input.hc_pbmc_pf_seuratobject_rds}" "{output.hc_pbmc_pf_b_seuratobject_rds}" &> "{log}"
    """

rule hc_pbmc_pf_myeloid_subsetting:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_pbmc_pf_myeloid_seuratobject_rds="output/subsets/hc_pbmc_pf_myeloid_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/hc_pbmc_pf_myeloid_subsetting.log",
  benchmark:
    "output/subsets/hc_pbmc_pf_myeloid_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/hc_pbmc_pf_myeloid_subsetting.R "{input.hc_pbmc_pf_seuratobject_rds}" "{output.hc_pbmc_pf_myeloid_seuratobject_rds}" &> "{log}"
    """

rule hc_pf_subsetting:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_pf_seuratobject_rds="output/subsets/hc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/hc_pf_subsetting.log",
  benchmark:
    "output/subsets/hc_pf_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript workflow/scripts/subsetting/hc_pf_subsetting.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{output.hc_pf_seuratobject_rds}" &> "{log}"
    """

## HC and CRC

rule hc_crc_pmp_pbmc_pf_subsetting:
  input:
    live_singlet_nonproliferating_seuratobject_rds="output/subsets/live_singlet_nonproliferating_SeuratObject.Rds",
  output:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_subsetting.log",
  benchmark:
    "output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_subsetting_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    txdonor="{txdonor}",
  shell:
    """
    Rscript workflow/scripts/subsetting/hc_crc_pmp_pbmc_pf_subsetting.R "{input.live_singlet_nonproliferating_seuratobject_rds}" "{params.txdonor}" "{output.hc_crc_pmp_pbmc_pf_seuratobject_rds}" &> "{log}"
    """

## CRC: PBMC, PF and TX


