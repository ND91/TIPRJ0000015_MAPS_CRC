
In the following analysis we will perform the single-cell RNA-sequencing analyses as described in XXX. To reproduce the analyses, you will need to do the following:
- Download the current snakemake pipeline using git (i.e. `git clone "XXX"`).
- Apply for acquiring the raw .fastq.gz and download the data from XXX (ID: XXX). Add the basepaths to the `config/samples/sample_files.xlsx`.
- Download Cellranger and add the basepath to `config/config.yml`. The preparation of the counts as described in the manuscript was done using Cellranger v7.0.0. 
- Download the Cellranger reference genome `refdata-gex-GRCh38-and-mm10-2020-A` and add the basepath to `config/config.yml`. 
- For some reason, certain Github R packages cannot be installed in a jobbed setting, you will need to install them manually by activating the environment and running `devtools::install_github`. This pertains the packages:
	- nichenetr
	- speckle
The way I approached this is to let the environment be created whereafter I activated it and manually installed the packages. Note that this must be redone anytime the environment is recreated for whatever reason.

Once you have completed the previous downloading the snakemake pipeline, create the conda environment to stage the subsequent steps.

```
conda env create -f environment.yml
``` 

Run the snakemake pipeline using the following command.

```
snakemake --cores X -p --use-conda
```

The `load=100` is meant as a way of controlling the number of instances to be run. For example, the `celltype_annotation` rule be greedy with memory. Left unchecked, the memory can run out.

The snakemake pipeline was run on my own machine using the following command.

```
snakemake --cores 20 -p --use-conda --resources mem_mb=189000 --rerun-triggers mtime
```