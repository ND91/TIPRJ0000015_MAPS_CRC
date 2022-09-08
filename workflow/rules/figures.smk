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