#!/usr/bin/env python3
# The goal of this script is to perform velocity analyses using scvelo on the PBMCs, PF, and TX monocytes and macrophages obtained from the same patients. 

import os
import argparse
import logging
import sys
import re
import anndata
import scvelo as scv

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)

def main(subset_anndata_h5ad, scvelo_anndata_h5ad):
    logging.info("Importing anndata")
    
    adata = anndata.read_h5ad(subset_anndata_h5ad)
    
    logging.info("Filtering genes")
    scv.pp.filter_genes(adata, min_shared_counts=20)
    
    logging.info("Normalizing cells")
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)
    
    logging.info("Calculating moments cells")
    scv.pp.moments(adata, n_pcs=43, n_neighbors=20)
    
    logging.info("Calculating velocity using dynamical model")
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.latent_time(adata)
    
    logging.info("Saving results in h5ad")
    adata.write_h5ad(scvelo_anndata_h5ad)
    
    return
    
if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--subset_anndata_h5ad", dest="subset_anndata_h5ad", required=True, nargs='+')
    argpar.add_argument("--scvelo_anndata_h5ad", dest="scvelo_anndata_h5ad", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.loom,
        parse_results.scvelo_loom,
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)

import scvelo as scv


