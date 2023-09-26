#!/usr/bin/env python3
# This script aims to subset the velocyto .loom file to monocytes and macrophages obtained from patients that provided PBMC, PF, and TX and store it as a anndata in h5ad.

import os
import shutil
import argparse
import logging
import sys
import re
import loompy
import numpy as np
import pandas as pd
import anndata

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)

def main(annotated_loom, subset_adata_h5ad):
    logging.info("Importing loom file")
    
    adata = anndata.read_loom(annotated_loom, validate = False)
    
    logging.info("Subsetting by TX donors")
    donor_tissue_pd = adata.obs[['Donor','Tissue']].drop_duplicates()
    txdonors_pd = donor_tissue_pd[donor_tissue_pd['Tissue'].isin(["TX"])] 
    adata_sub = adata[adata.obs['Donor'].isin(txdonors_pd['Donor'])]
    
    logging.info("Subsetting by monocytes and macrophages")
    adata_sub = adata_sub[adata_sub.obs['manual_l2'].isin(['Monocytes', 'Macrophages'])]
    
    logging.info("Saving anndata to h5ad")
    adata_sub.write_h5ad(subset_adata_h5ad)
    
    return
    
if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--subset_adata_h5ad", dest="subset_adata_h5ad", required=True)
    argpar.add_argument("--annotated_loom", dest="annotated_loom", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.annotated_loom,
        parse_results.subset_adata_h5ad,
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)
