# pbmc_pf_tx: 

## Differential analyses between tissues for all cells from all included patients.

rule pbmc_pf_tx_tissue_da:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/da/{pbmc_pf_tx_tissue_comparison}/{pbmc_pf_tx_tissue_comparison}_{relative_level}_pbmc_pf_tx_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx/da/{pbmc_pf_tx_tissue_comparison}/{pbmc_pf_tx_tissue_comparison}_{relative_level}_pbmc_pf_tx_propeller.log",
  benchmark:
    "output/analyses/pbmc_pf_tx/da/{pbmc_pf_tx_tissue_comparison}/{pbmc_pf_tx_tissue_comparison}_{relative_level}_pbmc_pf_tx_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    pbmc_pf_tx_tissue_comparison="{pbmc_pf_tx_tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_tissue.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.pbmc_pf_tx_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

rule pbmc_pf_tx_tissue_de:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/de/{pbmc_pf_tx_tissue_comparison}/{level}/{pbmc_pf_tx_tissue_comparison}_{level}_pbmc_pf_tx_degs_results.Rds",
    deseq2_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/de/{pbmc_pf_tx_tissue_comparison}/{level}/{pbmc_pf_tx_tissue_comparison}_{level}_pbmc_pf_tx_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx/de/{pbmc_pf_tx_tissue_comparison}/{level}/{pbmc_pf_tx_tissue_comparison}_{level}_pbmc_pf_tx_de.log",
  benchmark:
    "output/analyses/pbmc_pf_tx/de/{pbmc_pf_tx_tissue_comparison}/{level}/{pbmc_pf_tx_tissue_comparison}_{level}_pbmc_pf_tx_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    pbmc_pf_tx_tissue_comparison="{pbmc_pf_tx_tissue_comparison}",
    basepath="output/analyses/pbmc_pf_tx/de/{pbmc_pf_tx_tissue_comparison}/{level}/{pbmc_pf_tx_tissue_comparison}_{level}_pbmc_pf_tx"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_tissue.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.pbmc_pf_tx_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

## Differential analyses between peritoneal metastatic state for all cells from all included patients.

rule pbmc_pf_tx_metastasis_da:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_propeller.log",
  benchmark:
    "output/analyses/pbmc_pf_tx/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue="{tissue}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_metastasis.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

rule pbmc_pf_tx_metastasis_de:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_degs_results.Rds",
    deseq2_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_de.log",
  benchmark:
    "output/analyses/pbmc_pf_tx/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue="{tissue}",
    basepath="output/analyses/pbmc_pf_tx/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_metastasis.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

# pbmc_pf_tx_cleaned_cd45p: 

## Differential analyses between tissues for the cleaned CD45+ cells from all included patients.

rule pbmc_pf_tx_cleaned_cd45p_tissue_da:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_cd45p_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/da/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/da/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/da/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    pbmc_pf_tx_cleaned_cd45p_tissue_comparison="{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_tissue.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.pbmc_pf_tx_cleaned_cd45p_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

rule pbmc_pf_tx_cleaned_cd45p_tissue_de:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_cd45p_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/de/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{level}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{level}_pbmc_pf_tx_cleaned_cd45p_degs_results.Rds",
    deseq2_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/de/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{level}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{level}_pbmc_pf_tx_cleaned_cd45p_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/de/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{level}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{level}_pbmc_pf_tx_cleaned_cd45p_de.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/de/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{level}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{level}_pbmc_pf_tx_cleaned_cd45p_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    pbmc_pf_tx_cleaned_cd45p_tissue_comparison="{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}",
    basepath="output/analyses/pbmc_pf_tx_cleaned_cd45p/de/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}/{level}/{pbmc_pf_tx_cleaned_cd45p_tissue_comparison}_{level}_pbmc_pf_tx_cleaned_cd45p"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_tissue.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.pbmc_pf_tx_cleaned_cd45p_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

## Differential analyses between peritoneal metastatic state for the cleaned CD45+ cells from all included patients.

rule pbmc_pf_tx_cleaned_cd45p_metastasis_da:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_cd45p_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_pbmc_pf_tx_cleaned_cd45p_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue="{tissue}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_metastasis.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

rule pbmc_pf_tx_cleaned_cd45p_metastasis_de:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_cleaned_cd45p_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_cleaned_cd45p_degs_results.Rds",
    deseq2_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_cleaned_cd45p/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_cleaned_cd45p_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_cleaned_cd45p_de.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_cleaned_cd45p/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_cleaned_cd45p_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue="{tissue}",
    basepath="output/analyses/pbmc_pf_tx_cleaned_cd45p/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_pbmc_pf_tx_cleaned_cd45p"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_metastasis.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

# pbmc_pf_tx_paired: 

## Differential analyses between tissues for all cells from the paired patients.

rule pbmc_pf_tx_paired_tissue_da:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_paired_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_paired/da/{pbmc_pf_tx_paired_tissue_comparison}/{pbmc_pf_tx_paired_tissue_comparison}_{relative_level}_pbmc_pf_tx_paired_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_paired/da/{pbmc_pf_tx_paired_tissue_comparison}/{pbmc_pf_tx_paired_tissue_comparison}_{relative_level}_pbmc_pf_tx_paired_propeller.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_paired/da/{pbmc_pf_tx_paired_tissue_comparison}/{pbmc_pf_tx_paired_tissue_comparison}_{relative_level}_pbmc_pf_tx_paired_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    pbmc_pf_tx_paired_tissue_comparison="{pbmc_pf_tx_paired_tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_tissue.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.pbmc_pf_tx_paired_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

rule pbmc_pf_tx_paired_tissue_de:
  input:
    seuratobject_rds="output/subsets/pbmc_pf_tx_paired_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_paired/de/{pbmc_pf_tx_paired_tissue_comparison}/{level}/{pbmc_pf_tx_paired_tissue_comparison}_{level}_pbmc_pf_tx_paired_degs_results.Rds",
    deseq2_rds="/home/andrewliyim/ltytgat/MI-group/ayliyim/projects/TIPRJ0000015_MAPS/maps_crc/output/analyses/pbmc_pf_tx_paired/de/{pbmc_pf_tx_paired_tissue_comparison}/{level}/{pbmc_pf_tx_paired_tissue_comparison}_{level}_pbmc_pf_tx_paired_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/pbmc_pf_tx_paired/de/{pbmc_pf_tx_paired_tissue_comparison}/{level}/{pbmc_pf_tx_paired_tissue_comparison}_{level}_pbmc_pf_tx_paired_de.log",
  benchmark:
    "output/analyses/pbmc_pf_tx_paired/de/{pbmc_pf_tx_paired_tissue_comparison}/{level}/{pbmc_pf_tx_paired_tissue_comparison}_{level}_pbmc_pf_tx_paired_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    pbmc_pf_tx_paired_tissue_comparison="{pbmc_pf_tx_paired_tissue_comparison}",
    basepath="output/analyses/pbmc_pf_tx_paired/de/{pbmc_pf_tx_paired_tissue_comparison}/{level}/{pbmc_pf_tx_paired_tissue_comparison}_{level}_pbmc_pf_tx_paired"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_tissue.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.pbmc_pf_tx_paired_tissue_comparison}" "{input.functions_r}" &> "{log}"
    """
