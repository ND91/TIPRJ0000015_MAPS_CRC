In the following analysis we will perform the single-cell RNA-sequencing analyses as described in Saris et al. 2024. 

The `workflow` folder contains all the code necessary to prepare and analyze the data. Within the `workflow` folder, we have the following folders:
- `envs`: Conda/Mamba environments for running the scripts.
- `rules`: The Snakemake rules that orchestrate the different code to be run.
- `scripts`: The actual scripts that are being run.

Within the `scripts` folder, we have the following folders:
- `analyses`: Some scripts to show how we manually annotated our own data per lineage.
- `functions`: Folder containing helper functions.
- `manuscript`: Folder containing the scripts to prepare the figures for the manuscript and the review.
- `preparation`: Scripts for preparing our own data.
- `q1_pf_characterization`: Scripts belonging to part 1 where we interrogate the peritoneal immune system in healthy individuals.
- `q2_crc_vs_hc`: Scripts belonging to part 2 where we compare peritoneal immune system between patients with peritoneal metastasis from colorectal cancer origin with healthy individuals.
- `q3_pm_tx_characterization`: Scripts belonging to part 3 where we characterize the peritoneal metastatic immune system.
- `resources`: Scripts for preparing public data.
- `subsetting`: Script for subsetting the data.

To reproduce the analyses, you will need to do the following:
- Download the current snakemake pipeline using git (i.e. `git clone https://github.com/ND91/TIPRJ0000015_MAPS_CRC`).
- Apply for acquiring the raw .fastq.gz and download the data from EGA. PBMC and PF scRNA-seq and CITE-seq (PF) data from healthy controls can be found using accession IDs: EGAD50000000250, EGAD50000000251 and EGAD50000000252, respectively. PF and PM scRNA-seq data from PM-CRC patients can be found at accession IDs: EGAD50000000248 and EGAD50000000249, respectively. Add the basepaths to the `config/samples/sample_files.xlsx`.
- Acquire data from:
    - Liver Cell Atlas - Human (https://www.livercellatlas.org/download.php)
    - Gut Cell Atlas - Colon Immune Atlas (https://www.gutcellatlas.org/)
  Store it in `resources` 
- Download Cellranger and add the basepath to `config/config.yml`. The preparation of the counts as described in the manuscript was done using Cellranger v7.0.0. 
- Download the Cellranger reference genome `refdata-gex-GRCh38-and-mm10-2020-A` and add the basepath to `config/config.yml`. 
- For some reason, certain Github R packages cannot be installed in a jobbed setting, you will need to install them manually by activating the environment and running `devtools::install_github`. This pertains the packages:
	- speckle
The way I approached this is to let the environment be created whereafter I activated it and manually installed the packages. Note that this must be redone anytime the environment is recreated for whatever reason.

Important! The annotations of the macrophages in the code and the manuscript differ slightly:

Code: Manuscript
- Macrophages VCAN+: mono-CMs
- Macrophages C1Q+: C1Q+ CMs 
- Macrophages SPP1+: SPP1+ CMs 
- Macrophages VCAN+C1Q+: C1Q+ mono-CMs

Once you have completed the previous downloading the snakemake pipeline, create the conda/mamba environment (For running snakemake it is advisable to use `mamba`) to stage the subsequent steps.

```
mamba env create -f environment.yml
``` 

Run the snakemake pipeline using the following command.

```
snakemake --cores X -p --use-conda
```

The snakemake pipeline was run on my own machine using the following command.

```
snakemake --cores 20 -p --use-conda --resources mem_mb=189000 --rerun-triggers mtime
```

The aforementioned pipeline will generate more data than what was described in the manuscript (as is typical for scRNAseq analyses). After the necessary files are created, the scripts for generating the figures for the manuscript can be found in `workflow/scripts/manuscript` where each `.Rmd` file represents one of the main or supplementary figures. Note that some figures were generated in Graphpad Prism and that the final figures were assembled and edited for consistency purposes in Adobe Illustrator. The data generated for the reviewers can be found in `workflow/scripts/manuscript/review` where each file corresponds to a point raised by the reviewers.