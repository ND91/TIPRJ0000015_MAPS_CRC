# Differential abundance tissues: Compare all tissues, will be performed on both the unpaired and paired samples.
    
rule tissue_da:
  input:
    seuratobject_rds="output/subsets/{tissue_sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="output/analyses/{tissue_sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{tissue_sample_object}_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{tissue_sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{tissue_sample_object}_propeller.log",
  benchmark:
    "output/analyses/{tissue_sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{tissue_sample_object}_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue_comparison="{tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_tissue.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

# Differential abundance peritoneal metastasis: Compare all metastasis with non-metastasis, will be performed only on the unpaired samples.
    
rule metastasis_da:
  input:
    seuratobject_rds="output/subsets/{pm_sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="output/analyses/{pm_sample_object}/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_{pm_sample_object}_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{pm_sample_object}/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_{pm_sample_object}_propeller.log",
  benchmark:
    "output/analyses/{pm_sample_object}/da/PMpvPMn/PMpvPMn_{tissue}_{relative_level}_{pm_sample_object}_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue="{tissue}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_metastasis.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

# Differential expression tissues: Compare all tissues, will be performed on both the unpaired and paired samples.
    
rule tissue_de:
  input:
    seuratobject_rds="output/subsets/{tissue_sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="output/analyses/{tissue_sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{tissue_sample_object}_degs_results.Rds",
    deseq2_rds="output/analyses/{tissue_sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{tissue_sample_object}_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{tissue_sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{tissue_sample_object}_de.log",
  benchmark:
    "output/analyses/{tissue_sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{tissue_sample_object}_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue_comparison="{tissue_comparison}",
    basepath="output/analyses/{tissue_sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{tissue_sample_object}"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_tissue.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

# Differential expression peritoneal metastasis: Compare all metastasis with non-metastasis, will be performed only on the unpaired samples.
    
rule metastasis_de:
  input:
    seuratobject_rds="output/subsets/{pm_sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    degs_rds="output/analyses/{pm_sample_object}/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_{pm_sample_object}_degs_results.Rds",
    deseq2_rds="output/analyses/{pm_sample_object}/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_{pm_sample_object}_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{pm_sample_object}/de/{level}/PMpvPMn/PMpvPMn_{tissue}_{level}_{pm_sample_object}_de.log",
  benchmark:
    "output/analyses/{pm_sample_object}/de/{level}/PMpvPMn/PMpvPMn_{tissue}_{level}_{pm_sample_object}_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue="{tissue}",
    basepath="output/analyses/{pm_sample_object}/de/PMpvPMn/{level}/PMpvPMn_{tissue}_{level}_{pm_sample_object}"
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_metastasis.R "{input.seuratobject_rds}" "{output.degs_rds}" "{output.deseq2_rds}" "{params.basepath}" "{params.level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """
