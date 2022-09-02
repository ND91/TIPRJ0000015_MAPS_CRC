rule fig_cells_per_organism:
  input:
    counts=expand("output/cellranger/{runid}/outs/raw_feature_bc_matrix", runid = runs),
    sample_metadata=config['sample_metadata'],
  output:
    cells_per_organism_raw_rds="output/figures/cells_per_organism/cells_per_organism_raw_df.Rds",
    cells_per_organism_raw_csv="output/figures/cells_per_organism/cells_per_organism_raw_df.csv",
    fig_cells_per_organism_raw_pdf="output/figures/cells_per_organism/cells_per_organism_raw.pdf",
    cells_per_organism_processed_rds="output/figures/cells_per_organism/cells_per_organism_processed_df.Rds",
    cells_per_organism_processed_csv="output/figures/cells_per_organism/cells_per_organism_processed_df.csv",
    fig_cells_per_organism_processed_pdf="output/figures/cells_per_organism/cells_per_organism_processed.pdf",
  threads: 
    1
  conda:
    "../envs/seurat.yaml"
  log:
    "output/figures/cells_per_organism.log",
  benchmark:
    "output/figures/cells_per_organism_benchmark.txt"
  resources:
    mem_mb=150000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_cells_per_organism.R {input.counts} {input.sample_metadata} {output.cells_per_organism_raw_rds} {output.cells_per_organism_raw_csv} {output.fig_cells_per_organism_raw_pdf} {output.cells_per_organism_processed_rds} {output.cells_per_organism_processed_csv} {output.fig_cells_per_organism_processed_pdf} &>{log}
    """

rule cells_per_sampleID

rule cells_per_runID

rule umap_spleen_pbmc:
  input:
    spleen_pbmc_seuratobject_rds="output/spleen_pbmc/spleen_pbmc_SeuratObject.Rds",
  output:
    umap_spleen_pbmc_pdf="output/figures/umap_spleen_pbmc.pdf",
  threads: 
    1
  conda:
    "../envs/seurat.yaml",
  log:
    "output/figures/umap_hs_cd45p_live_singlet_spleen_pbmc.log",
  benchmark:
    "output/figures/umap_hs_cd45p_live_singlet_spleen_pbmc_benchmark.txt",
  resources:
    mem_mb=60000,
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_hs_cd45p_live_singlet_wrefpbmc_umap_pbmc_spleen.R {input.hs_cd45p_spleen_pbmc_seuratobject_rds} {output.umap_hs_cd45p_live_singlet_spleen_pbmc_pdf} &> {log}
    """

rule hs_cd45p_live_singlet_wrefpbmc_umap_pbmc

rule hs_cd45p_live_singlet_wrefpbmc_umap_spleen

rule_hs_cd45p_live_singlet_wrefpbmc_abundance_boxplot

rule_hs_cd45p_live_singlet_wrefpbmc_abundance_heatmap

rule_hs_cd45p_live_singlet_wrefpbmc_expression_heatmap

