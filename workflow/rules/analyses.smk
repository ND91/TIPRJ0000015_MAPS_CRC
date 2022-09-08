# Differential abundance tissues
    
rule tissue_da:
  input:
    seuratobject_rds="output/subsets/{sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="output/analyses/{sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{sample_object}_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{sample_object}_propeller.log",
  benchmark:
    "output/analyses/{sample_object}/da/{tissue_comparison}/{tissue_comparison}_{relative_level}_{sample_object}_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue_comparison="{tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_tissue.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

# Differential abundance peritoneal metastasis
    
rule metastasis_da:
  input:
    seuratobject_rds="output/subsets/{sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    propeller_rds="output/analyses/{sample_object}/da/pmpvpmn/pmpvpmn_{tissue}_{relative_level}_{sample_object}_propeller_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{sample_object}/da/pmpvpmn/pmpvpmn_{tissue}_{relative_level}_{sample_object}_propeller.log",
  benchmark:
    "output/analyses/{sample_object}/da/pmpvpmn/pmpvpmn_{tissue}_{relative_level}_{sample_object}_propeller_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    relative_level="{relative_level}",
    tissue="{tissue}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/da_analyses/da_metastasis.R "{input.seuratobject_rds}" "{output.propeller_rds}" "{params.relative_level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

# Differential expression tissues
    
rule tissue_de:
  input:
    seuratobject_rds="output/subsets/{sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    deseq2_rds="output/analyses/{sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{sample_object}_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{sample_object}_de.log",
  benchmark:
    "output/analyses/{sample_object}/de/{tissue_comparison}/{level}/{tissue_comparison}_{level}_{sample_object}_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue_comparison="{tissue_comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_tissue.R "{input.seuratobject_rds}" "{output.deseq2_rds}" "{params.level}" "{params.tissue_comparison}" "{input.functions_r}" &> "{log}"
    """

# Differential expression peritoneal metastasis
    
rule metastasis_de:
  input:
    seuratobject_rds="output/subsets/{sample_object}_SeuratObject.Rds",
    functions_r="workflow/scripts/functions.R",
  output:
    deseq2_rds="output/analyses/{sample_object}/de/pmpvpmn/{level}/pmpvpmn_{tissue}_{level}_{sample_object}_deseq2_results.Rds",
  threads: 
    1
  conda:
    "../envs/deseq2.yaml",
  log:
    "output/analyses/{sample_object}/de/pmpvpmn/{level}/pmpvpmn_{tissue}_{level}_{sample_object}_de.log",
  benchmark:
    "output/analyses/{sample_object}/de/pmpvpmn/{level}/pmpvpmn_{tissue}_{level}_{sample_object}_de_benchmark.txt",
  resources:
    mem_mb=60000,
  params:
    level="{level}",
    tissue="{tissue}",
  shell:
    """
    Rscript --vanilla workflow/scripts/analyses/de_analyses/pseudobulk_de_metastasis.R "{input.seuratobject_rds}" "{output.deseq2_rds}" "{params.level}" "{params.tissue}" "{input.functions_r}" &> "{log}"
    """

