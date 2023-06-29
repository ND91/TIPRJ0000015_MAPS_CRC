# HC

rule fig_hc_umap_pbmc_pf_coltissue:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_umap_pbmc_pf_coltissue_pdf="output/figures/hc/hc_umap_pbmc_pf_coltissue.pdf",
    hc_umap_pbmc_pf_coltissue_splittissue_pdf="output/figures/hc/hc_umap_pbmc_pf_coltissue_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_umap_pbmc_pf_coltissue.log",
  benchmark:
    "output/figures/hc/hc_umap_pbmc_pf_coltissue_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_umap_pbmc_pf_coltissue.R "{input.hc_pbmc_pf_seuratobject_rds}" "{output.hc_umap_pbmc_pf_coltissue_pdf}" "{output.hc_umap_pbmc_pf_coltissue_splittissue_pdf}" &> "{log}"
    """

rule fig_hc_umap_pbmc_pf_colcelltype_splittissue:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_umap_pbmc_pf_colcelltype_splittissue_pdf="output/figures/hc/hc_umap_pbmc_pf_col{level}_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_umap_pbmc_pf_col{level}_splittissue.log",
  benchmark:
    "output/figures/hc/hc_umap_pbmc_pf_col{level}_splittissue_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_umap_pbmc_pf_colcelltype_splittissue.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_umap_pbmc_pf_colcelltype_splittissue_pdf}" &> "{log}"
    """

rule fig_hc_umap_pf_colcelltype:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_gexumap_pf_coltissue_pdf="output/figures/hc/hc_umap_gex_pf_col{level}.pdf",
    hc_citeumap_pf_coltissue_pdf="output/figures/hc/hc_umap_cite_pf_col{level}.pdf",
    hc_wnnumap_pf_coltissue_pdf="output/figures/hc/hc_umap_wnn_pf_col{level}.pdf",
    hc_totalumap_pf_coltissue_pdf="output/figures/hc/hc_umap_total_pf_col{level}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_umap_pf_col{level}.log",
  benchmark:
    "output/figures/hc/hc_umap_pf_col{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_umap_pf_colcelltype.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_gexumap_pf_coltissue_pdf}" "{output.hc_citeumap_pf_coltissue_pdf}" "{output.hc_wnnumap_pf_coltissue_pdf}" "{output.hc_totalumap_pf_coltissue_pdf}" &> "{log}"
    """

rule fig_hc_dotplot_markerexpression_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_dotplot_markerexpression_pbmc_pf_pdf="output/figures/hc/hc_dotplot_markerexpression_pbmc_pf_{level}.pdf",
    hc_dotplot_markerexpression_pbmc_pf_splittissue_pdf="output/figures/hc/hc_dotplot_markerexpression_pbmc_pf_{level}_splittissue.pdf",
    hc_dotplot_markerexpression_pf_gex_pdf="output/figures/hc/hc_dotplot_markerexpression_pf_{level}_gex.pdf",
    hc_dotplot_markerexpression_pf_pex_pdf="output/figures/hc/hc_dotplot_markerexpression_pf_{level}_pex.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_dotplot_markerexpression_pbmc_pf_{level}.log",
  benchmark:
    "output/figures/hc/hc_dotplot_markerexpression_pbmc_pf_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_dotplot_markerexpression_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_dotplot_markerexpression_pbmc_pf_pdf}" "{output.hc_dotplot_markerexpression_pbmc_pf_splittissue_pdf}" "{output.hc_dotplot_markerexpression_pf_gex_pdf}" "{output.hc_dotplot_markerexpression_pf_pex_pdf}" &> "{log}"
    """

rule fig_hc_heatmap_markerexpression_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_heatmap_markerexpression_pbmc_pf_pdf="output/figures/hc/hc_heatmap_markerexpression_pbmc_pf_{level}.pdf",
    hc_heatmap_markerexpression_pbmc_pf_splittissue_pdf="output/figures/hc/hc_heatmap_markerexpression_pbmc_pf_{level}_splittissue.pdf",
    hc_heatmap_markerexpression_pf_gex_pdf="output/figures/hc/hc_heatmap_markerexpression_pf_{level}_gex.pdf",
    hc_heatmap_markerexpression_pf_pex_pdf="output/figures/hc/hc_heatmap_markerexpression_pf_{level}_pex.pdf",
    hc_heatmap_markerexpression_pf_gex_pex_pdf="output/figures/hc/hc_heatmap_markerexpression_pf_{level}_gex_pex.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/hc/hc_heatmap_markerexpression_pbmc_pf_{level}.log",
  benchmark:
    "output/figures/hc/hc_heatmap_markerexpression_pbmc_pf_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_heatmap_markerexpression_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_heatmap_markerexpression_pbmc_pf_pdf}" "{output.hc_heatmap_markerexpression_pbmc_pf_splittissue_pdf}" "{output.hc_heatmap_markerexpression_pf_gex_pdf}" "{output.hc_heatmap_markerexpression_pf_pex_pdf}" "{output.hc_heatmap_markerexpression_pf_gex_pex_pdf}" &> "{log}"
    """

rule fig_hc_violinplot_pbmc_pf:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_violinplot_pbmc_pf_splittissue_pdf="output/figures/hc/hc_violinplot_pbmc_pf_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_violinplot_pbmc_pf.log",
  benchmark:
    "output/figures/hc/hc_violinplot_pbmc_pf_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_violinplot_pbmc_pf.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_violinplot_pbmc_pf_splittissue_pdf}" &> "{log}"
    """

rule fig_hc_pfvpbmc_de:
  input:
    deseq2_list_rds="output/analyses/hc_pfvpbmc/hc_pfvpbmc_{level}_deseq2_list.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    arrowplot_top5sig_pdf="output/figures/hc/hc_pfvpbmc_{level}_arrowplot_top5sig.pdf",
    arrowplot_resident_pdf="output/figures/hc/hc_pfvpbmc_{level}_arrowplot_resident.pdf",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/figures/hc/hc_pfvpbmc_{level}_de.log",
  benchmark:
    "output/figures/hc/hc_pfvpbmc_{level}_de_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_pfvpbmc_de.R "{input.deseq2_list_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.arrowplot_top5sig_pdf}" "{output.arrowplot_resident_pdf}" &> "{log}"
    """

rule fig_hc_boxplot_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_boxplot_pbmc_abundance_pdf="output/figures/hc/hc_boxplot_pbmc_abundance_{comparison}.pdf",
    hc_boxplot_pbmc_abundance_patanno_pdf="output/figures/hc/hc_boxplot_pbmc_abundance_patanno_{comparison}.pdf",
    hc_boxplot_pf_abundance_pdf="output/figures/hc/hc_boxplot_pf_abundance_{comparison}.pdf",
    hc_boxplot_pf_abundance_patanno_pdf="output/figures/hc/hc_boxplot_pf_abundance_patanno_{comparison}.pdf",
    hc_boxplot_pbmc_pf_abundance_pdf="output/figures/hc/hc_boxplot_pbmc_pf_abundance_{comparison}.pdf",
    hc_boxplot_pbmc_pf_abundance_patanno_pdf="output/figures/hc/hc_boxplot_pbmc_pf_abundance_patanno_{comparison}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_boxplot_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/figures/hc/hc_boxplot_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_boxplot_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.hc_boxplot_pbmc_abundance_pdf}" "{output.hc_boxplot_pbmc_abundance_patanno_pdf}" "{output.hc_boxplot_pf_abundance_pdf}" "{output.hc_boxplot_pf_abundance_patanno_pdf}" "{output.hc_boxplot_pbmc_pf_abundance_pdf}" "{output.hc_boxplot_pbmc_pf_abundance_patanno_pdf}" &> "{log}"
    """

rule fig_hc_boxplot_pbmc_pf_liver_colon_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    hc_liver_seuratobject_rds="resources/liver/liver_seuratobject.rds",
    hc_colon_seuratobject_rds="resources/colon/colon_seuratobject.rds",
    marker_genes_xlsx=config['celltype_markers'],
  output:
    hc_boxplot_pbmc_pf_colon_liver_abundance_pdf="output/figures/hc/hc_boxplot_pbmc_pf_colon_liver_abundance.pdf",
    hc_boxplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf="output/figures/hc/hc_boxplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_boxplot_pbmc_pf_liver_colon_abundance.log",
  benchmark:
    "output/figures/hc/hc_boxplot_pbmc_pf_liver_colon_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_boxplot_pbmc_pf_liver_colon_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.hc_liver_seuratobject_rds}" "{input.hc_colon_seuratobject_rds}" "{input.marker_genes_xlsx}" "{output.hc_boxplot_pbmc_pf_colon_liver_abundance_pdf}" "{output.hc_boxplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf}" &> "{log}"
    """

rule fig_hc_stackedbarplot_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_stackedbarplot_pbmc_abundance_pdf="output/figures/hc/hc_stackedbarplot_pbmc_abundance_{level}.pdf",
    hc_stackedbarplot_pf_abundance_pdf="output/figures/hc/hc_stackedbarplot_pf_abundance_{level}.pdf",
    hc_stackedbarplot_pbmc_pf_abundance_pdf="output/figures/hc/hc_stackedbarplot_pbmc_pf_abundance_{level}.pdf",
    hc_stackedbarplot_pbmc_abundance_myeloid_pdf="output/figures/hc/hc_stackedbarplot_pbmc_abundance_{level}_myeloid.pdf",
    hc_stackedbarplot_pf_abundance_myeloid_pdf="output/figures/hc/hc_stackedbarplot_pf_abundance_{level}_myeloid.pdf",
    hc_stackedbarplot_pbmc_pf_abundance_myeloid_pdf="output/figures/hc/hc_stackedbarplot_pbmc_pf_abundance_{level}_myeloid.pdf",
    hc_stackedbarplot_pbmc_abundance_t_pdf="output/figures/hc/hc_stackedbarplot_pbmc_abundance_{level}_t.pdf",
    hc_stackedbarplot_pf_abundance_t_pdf="output/figures/hc/hc_stackedbarplot_pf_abundance_{level}_t.pdf",
    hc_stackedbarplot_pbmc_pf_abundance_t_pdf="output/figures/hc/hc_stackedbarplot_pbmc_pf_abundance_{level}_t.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_stackedbarplot_pbmc_pf_abundance_{level}.log",
  benchmark:
    "output/figures/hc/hc_stackedbarplot_pbmc_pf_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_stackedbarplot_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.level}" "{output.hc_stackedbarplot_pbmc_abundance_pdf}" "{output.hc_stackedbarplot_pf_abundance_pdf}" "{output.hc_stackedbarplot_pbmc_pf_abundance_pdf}" "{output.hc_stackedbarplot_pbmc_abundance_myeloid_pdf}" "{output.hc_stackedbarplot_pf_abundance_myeloid_pdf}" "{output.hc_stackedbarplot_pbmc_pf_abundance_myeloid_pdf}" "{output.hc_stackedbarplot_pbmc_abundance_t_pdf}" "{output.hc_stackedbarplot_pf_abundance_t_pdf}" "{output.hc_stackedbarplot_pbmc_pf_abundance_t_pdf}" &> "{log}"
    """

rule fig_hc_stackedbarplot_pbmc_pf_liver_colon_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
    hc_colon_seuratobject_rds="resources/colon/colon_seuratobject.rds",
    hc_liver_seuratobject_rds="resources/liver/liver_seuratobject.rds",
    marker_genes_xlsx=config['celltype_markers'],
  output:
    hc_stackedbarplot_pbmc_pf_colon_liver_abundance_pdf="output/figures/hc/hc_stackedbarplot_pbmc_pf_colon_liver_abundance.pdf",
    hc_stackedbarplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf="output/figures/hc/hc_stackedbarplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_stackedbarplot_pbmc_pf_liver_colon_abundance.log",
  benchmark:
    "output/figures/hc/hc_stackedbarplot_pbmc_pf_liver_colon_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_stackedbarplot_pbmc_pf_liver_colon_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{input.hc_colon_seuratobject_rds}" "{input.hc_liver_seuratobject_rds}" "{input.marker_genes_xlsx}" "{output.hc_stackedbarplot_pbmc_pf_colon_liver_abundance_pdf}" "{output.hc_stackedbarplot_pbmc_pf_caecum_transverse_sigmoid_liver_abundance_pdf}" &> "{log}"
    """

rule fig_hc_heatmap_pbmc_pf_abundance:
  input:
    hc_pbmc_pf_seuratobject_rds="output/subsets/hc_pbmc_pf_SeuratObject.Rds",
  output:
    hc_heatmapplot_pbmc_pf_abundance_scaled_pdf="output/figures/hc/hc_heatmap_pbmc_pf_abundance_{comparison}_scaled.pdf",
    hc_heatmapplot_pbmc_pf_abundance_unscaled_pdf="output/figures/hc/hc_heatmap_pbmc_pf_abundance_{comparison}_unscaled.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/hc/hc_heatmap_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/figures/hc/hc_heatmap_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_heatmap_pbmc_pf_abundance.R "{input.hc_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.hc_heatmapplot_pbmc_pf_abundance_scaled_pdf}" "{output.hc_heatmapplot_pbmc_pf_abundance_unscaled_pdf}" &> "{log}"
    """

rule fig_hc_mnp_markerexpression_pf:
  input:
    hc_pf_mnp_seuratobject_rds="output/subsets/hc_pf_mnp_SeuratObject.Rds",
  output:
    hc_violinplot_mnp_markerexpression_pf_pdf="output/figures/hc/hc_violinplot_mnp_markerexpression_pf.pdf",
    hc_boxplot_mnp_markerexpression_pf_pdf="output/figures/hc/hc_boxplot_mnp_markerexpression_pf.pdf",
    hc_dotplot_mnp_markerexpression_pf_pdf="output/figures/hc/hc_dotplot_mnp_markerexpression_pf.pdf",
    hc_umap_mnp_markerexpression_pf_pdf="output/figures/hc/hc_umap_mnp_markerexpression_pf.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc/hc_mnp_markerexpression_pf.log",
  benchmark:
    "output/figures/hc/hc_mnp_markerexpression_pf_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_mnp_markerexpression_pf.R "{input.hc_pf_mnp_seuratobject_rds}" "{output.hc_violinplot_mnp_markerexpression_pf_pdf}" "{output.hc_boxplot_mnp_markerexpression_pf_pdf}" "{output.hc_dotplot_mnp_markerexpression_pf_pdf}" "{output.hc_umap_mnp_markerexpression_pf_pdf}" &> "{log}"
    """

# HC & CRC

rule fig_hc_crc_pmp_umap_pf_colcelltype_splitgroup:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    hc_crc_pmp_umap_pf_colcelltype_splitgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_col{level}_umap_pf_colcelltype_splitgroup_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_umap_pf_colcelltype_splitgroup.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.hc_crc_pmp_umap_pf_colcelltype_splitgroup_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_umap_pf_colgroup:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_umap_pf_colgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_umap_pf_colgroup.pdf",
    hc_crc_pmp_umap_pf_colgroup_splitgroup_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_umap_pf_colgroup_splitgroup.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/fig_hc_crc_pmp_{txdonor}_umap_pf_colgroup.log",
  benchmark:
    "output/figures/hc_crc/fig_hc_crc_pmp_{txdonor}_umap_pf_colgroup_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_umap_pf_colgroup.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{output.hc_crc_pmp_umap_pf_colgroup_pdf}" "{output.hc_crc_pmp_umap_pf_colgroup_splitgroup_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_boxplot_pbmc_pf_abundance:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_boxplot_pbmc_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_abundance_{comparison}.pdf",
    hc_crc_pmp_boxplot_pf_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_abundance_{comparison}.pdf",
    hc_crc_pmp_boxplot_pbmc_pf_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pbmc_pf_abundance_{comparison}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_boxplot_pbmc_pf_abundance.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{params.comparison}" "{output.hc_crc_pmp_boxplot_pbmc_abundance_pdf}" "{output.hc_crc_pmp_boxplot_pf_abundance_pdf}" "{output.hc_crc_pmp_boxplot_pbmc_pf_abundance_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance:
  input:
    hc_crc_pmp_pbmc_pf_seuratobject_rds="output/subsets/hc_crc_pmp_{txdonor}_pbmc_pf_SeuratObject.Rds",
  output:
    hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_boxplot_pf_macrophage_l4rl2_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance.R "{input.hc_crc_pmp_pbmc_pf_seuratobject_rds}" "{output.hc_crc_pmp_boxplot_pf_macrophage_l4rl2_abundance_pdf}" &> "{log}"
    """

rule fig_hc_crc_pmp_volcanoplot_pf_macrophage_l3_crcpmpvhc:
  input:
    deseq2_list_rds="output/analyses/hc_crc_pmp_{txdonor}_pf_macrophages_l3_crcpmpvhc_deseq2_list.Rds",
  output:
    volcanoplot_pdf="output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_crcpmpvhc.pdf",
  threads: 
    1
  conda:
    "../envs/r-deseq2.yaml",
  log:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_macrophage_l3_crcpmpvhc.log",
  benchmark:
    "output/figures/hc_crc/hc_crc_pmp_{txdonor}_volcanoplot_pf_macrophage_l3_crcpmpvhc_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hc_crc_pmp_volcanoplot_pf_macrophages_l3_crcpmpvhc_de.R "{input.deseq2_list_rds}" "{output.volcanoplot_pdf}" &> "{log}"
    """

# CRC

rule fig_crc_patient_metadata:
  input:
    patient_metadata_xlsx=config["patient_metadata"],
  output:
    crc_heatmap_patient_metadata_pdf="output/figures/crc/crc_heatmap_patient_metadata.pdf",
    crc_heatmap_patient_metadata_nocytof_pdf="output/figures/crc/crc_heatmap_patient_metadata_nocytof.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/crc/crc_patient_metadata.log",
  benchmark:
    "output/figures/crc/crc_patient_metadata_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_heatmap_patient_metadata.R "{input.patient_metadata_xlsx}" "{output.crc_heatmap_patient_metadata_pdf}" "{output.crc_heatmap_patient_metadata_nocytof_pdf}" &> "{log}"
    """ 

rule fig_crc_umap_pbmc_pf_tx_colcelltype_splittissue:
  input:
    crc_pbmc_pf_tx_paired_seuratobject_rds="output/subsets/crc_pbmc_pf_tx_paired_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_umap_pbmc_pf_tx_paired_colcelltype_splittissue_pdf="output/figures/crc/crc_umap_pbmc_pf_tx_paired_col{level}_splittissue.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_umap_pbmc_pf_tx_col{level}_splittissue.log",
  benchmark:
    "output/figures/crc/crc_umap_pbmc_pf_tx_col{level}_splittissue_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_umap_pbmc_pf_tx_colcelltype_splittissue.R "{input.crc_pbmc_pf_tx_paired_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_umap_pbmc_pf_tx_paired_colcelltype_splittissue_pdf}" &> "{log}"
    """

rule fig_crc_boxplot_tx_abundance:
  input:
    crc_tx_seuratobject_rds="output/subsets/crc_tx_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_boxplot_tx_abundance_top5_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_top5.pdf",
    crc_boxplot_tx_abundance_top10_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_top10.pdf",
    crc_boxplot_tx_abundance_all_pdf="output/figures/crc/crc_boxplot_tx_abundance_{level}_all.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_tx_abundance_{level}.log",
  benchmark:
    "output/figures/crc/crc_boxplot_tx_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_boxplot_tx_abundance.R "{input.crc_tx_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_boxplot_tx_abundance_top5_pdf}" "{output.crc_boxplot_tx_abundance_top10_pdf}" "{output.crc_boxplot_tx_abundance_all_pdf}" &> "{log}"
    """

rule fig_crc_boxplot_myeloid_tx_abundance:
  input:
    crc_myeloid_tx_seuratobject_rds="output/subsets/crc_tx_myeloid_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_boxplot_myeloid_tx_abundance_top5_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_top5.pdf",
    crc_boxplot_myeloid_tx_abundance_top10_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_top10.pdf",
    crc_boxplot_myeloid_tx_abundance_all_pdf="output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_all.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}.log",
  benchmark:
    "output/figures/crc/crc_boxplot_myeloid_tx_abundance_{level}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_boxplot_myeloid_tx_abundance.R "{input.crc_myeloid_tx_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{params.level}" "{output.crc_boxplot_myeloid_tx_abundance_top5_pdf}" "{output.crc_boxplot_myeloid_tx_abundance_top10_pdf}" "{output.crc_boxplot_myeloid_tx_abundance_all_pdf}" &> "{log}"
    """

rule fig_crc_scatterplot_pf_tx_macrophage_l4rl2_pci:
  input:
    crc_pbmc_pf_tx_paired_seuratobject_rds="output/subsets/crc_pbmc_pf_tx_paired_SeuratObject.Rds",
    celltype_markers_xlsx=config['celltype_markers'],
  output:
    crc_scatterplot_pf_tx_macrophage_l4rl2_pci_pdf="output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci.log",
  benchmark:
    "output/figures/crc/crc_scatterplot_pf_tx_macrophage_l4rl2_pci_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    level="{level}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_scatterplot_macrophage_l4_pci.R "{input.crc_pbmc_pf_tx_paired_seuratobject_rds}" "{input.celltype_markers_xlsx}" "{output.crc_scatterplot_pf_tx_macrophage_l4rl2_pci_pdf}" &> "{log}"
    """

rule fig_crc_t_clonotype_abundance:
  input:
    crc_pbmc_pf_tx_t_paired_cdr3_airr_csv="output/trust4/crc_pbmc_pf_tx_t_paired_airr.csv",
  output:
    crc_upset_t_clonotype_abundance_pdf="output/figures/crc/crc_upset_t_clonotype_abundance_{patient}.pdf",
  threads: 
    1
  conda:
    "../envs/r-complexheatmap.yaml",
  log:
    "output/figures/crc/crc_t_clonotype_abundance_{patient}.log",
  benchmark:
    "output/figures/crc/crc_t_clonotype_abundance_{patient}_benchmark.txt",
  resources:
    mem_mb=16000,
  params:
    patient="{patient}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_t_clonotype_abundance.R "{input.crc_pbmc_pf_tx_t_paired_cdr3_airr_csv}" "{params.patient}" "{output.crc_upset_t_clonotype_abundance_pdf}" &> "{log}"
    """

rule fig_crc_boxplot_t_pairwise_clonotype_abundance:
  input:
    crc_pbmc_pf_tx_t_paired_cdr3_airr_csv="output/trust4/crc_pbmc_pf_tx_t_paired_airr.csv",
  output:
    crc_upset_t_clonotype_abundance_pdf="output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance.pdf",
    crc_upset_t_clonotype_abundance_patanno_pdf="output/figures/crc/crc_boxplot_t_pairwise_clonotype_patanno_abundance.pdf",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance.log",
  benchmark:
    "output/figures/crc/crc_boxplot_t_pairwise_clonotype_abundance_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_crc_boxplot_t_pairwise_clonotype_abundance.R "{input.crc_pbmc_pf_tx_t_paired_cdr3_airr_csv}" "{output.crc_upset_t_clonotype_abundance_pdf}" "{output.crc_upset_t_clonotype_abundance_patanno_pdf}" &> "{log}"
    """










rule fig_boxplot_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/subsets/tx_SeuratObject.Rds",
  output:
    boxplot_tx_l3rl2_df_csv="output/figures/fig_boxplot_tx_l3rl2.csv",
    boxplot_tx_l3rl2_svg="output/figures/fig_boxplot_tx_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_tx_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_tx_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_tx_l3rl2.R {input.tx_seuratobject_rds} {output.boxplot_tx_l3rl2_df_csv} {output.boxplot_tx_l3rl2_svg} &> {log}
    """

rule fig_tsne_pbmc_pf_tx_tsne:
  input:
    pbmc_pf_tx_paired_seuratobject_rds="output/subsets/pbmc_pf_tx_paired_SeuratObject.Rds",
  output:
    tsne_pbmc_pf_tx_df_csv="output/figures/tsne_pbmc_pf_tx.csv",
    tsne_pbmc_pf_tx_fig1_svg="output/figures/fig_tsne_pbmc_pf_tx_fig1.svg",
    tsne_pbmc_pf_tx_fig2_svg="output/figures/fig_tsne_pbmc_pf_tx_fig2.svg",
    tsne_pbmc_pf_tx_fig3_svg="output/figures/fig_tsne_pbmc_pf_tx_fig3.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_pbmc_pf_tx_tsne.log",
  benchmark:
    "output/figures/fig_pbmc_pf_tx_tsne_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_tsne_pbmc_pf_tx.R {input.pbmc_pf_tx_paired_seuratobject_rds} {output.tsne_pbmc_pf_tx_df_csv} {output.tsne_pbmc_pf_tx_fig1_svg} {output.tsne_pbmc_pf_tx_fig2_svg} {output.tsne_pbmc_pf_tx_fig3_svg} &> {log}
    """

rule fig_boxplot_pf_tx_l3rl2:
  input:
    tx_seuratobject_rds="output/subsets/tx_SeuratObject.Rds",
  output:
    boxplot_pf_tx_l3rl2_df_csv="output/figures/fig_boxplot_pf_tx_l3rl2.csv",
    boxplot_pf_tx_l3rl2_svg="output/figures/fig_boxplot_pf_tx_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_pf_tx_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_pf_tx_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_pf_tx_l3rl2.R "{input.tx_seuratobject_rds}" "{output.boxplot_pf_tx_l3rl2_df_csv}" "{output.boxplot_pf_tx_l3rl2_svg}" &> "{log}"
    """

rule fig_boxplot_pbmc_pf_l3rl2:
  input:
    pbmc_pf_tx_seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
  output:
    boxplot_pbmc_pf_l3rl2_df_csv="output/figures/fig_boxplot_pbmc_pf_l3rl2.csv",
    boxplot_pbmc_pf_l3rl2_svg="output/figures/fig_boxplot_pbmc_pf_l3rl2.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_boxplot_pbmc_pf_l3rl2.log",
  benchmark:
    "output/figures/fig_boxplot_pbmc_pf_l3rl2_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_boxplot_pbmc_pf_l3rl2.R "{input.pbmc_pf_tx_seuratobject_rds}" "{output.boxplot_pbmc_pf_l3rl2_df_csv}" "{output.boxplot_pbmc_pf_l3rl2_svg}" &> "{log}"
    """

rule fig_heatmap_cd4t_markers:
  input:
    pbmc_pf_tx_seuratobject_rds="output/subsets/pbmc_pf_tx_SeuratObject.Rds",
    markers_csv=config["markers"],
  output:
    heatmap_cd4t_markers_svg="output/figures/fig_heatmap_cd4t_markers.svg",
  threads: 
    1
  conda:
    "../envs/r.yaml",
  log:
    "output/figures/fig_heatmap_cd4t_markers.log",
  benchmark:
    "output/figures/fig_heatmap_cd4t_markers_benchmark.txt",
  resources:
    mem_mb=16000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_heatmap_markers_cd4t.R {input.pbmc_pf_tx_seuratobject_rds} {input.markers_csv} {output.heatmap_cd4t_markers_svg} &> {log}
    """
