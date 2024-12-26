#!/usr/bin/env python3
# This script aims to perform scvelo analyses.

import os
import shutil
import argparse
import logging
import sys
import re
import loompy
import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)

def main(cellmetadata_csv, counts_mtx, features_csv, umap_csv, pca_csv, velocyto_curated_loom, scvelo_anndata_h5ad, scvelo_anndata_cellmetadata_csv):
    logging.info("Importing flat files")
    cellmetadata = pd.read_csv(cellmetadata_csv)
    counts = io.mmread(counts_mtx)
    with open(features_csv, 'r') as f:
        features = f.read().splitlines()
    pca = pd.read_csv(pca_csv)
    umap = pd.read_csv(umap_csv)

    logging.info("Creating AnnData")
    adata = anndata.AnnData(
        X=counts.transpose().tocsr()
    )
    adata.obs = cellmetadata
    adata.obs.index = adata.obs['CellID']
    adata.var.index = features
    
    umap.index = pca.index = adata.obs.index
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = umap.to_numpy()
    
    logging.info("Importing velocyto loom file")
    velocyto_curated_anndata = anndata.read_loom(velocyto_curated_loom, validate = False)
    
    logging.info("Creating subsetted anndata file")
    #velocyto_curated_subset_anndata = velocyto_curated_anndata[velocyto_curated_anndata.obs.index.isin(adata.obs.index)]
    velocyto_curated_subset_anndata = velocyto_curated_anndata[adata.obs.index] # I can see this going wrong if there are cells in the subsetted adata and not in the velocyto_curated.
    velocyto_curated_subset_anndata.obsm['X_pca'] = adata.obsm['X_pca']
    velocyto_curated_subset_anndata.obsm['X_umap'] = adata.obsm['X_umap']
    velocyto_curated_subset_anndata.obs = adata.obs
  
    logging.info("Performing scvelo analyses")
    scv.pp.filter_and_normalize(velocyto_curated_subset_anndata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(velocyto_curated_subset_anndata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(velocyto_curated_subset_anndata)
    scv.tl.velocity_graph(velocyto_curated_subset_anndata)
    scv.tl.velocity_pseudotime(velocyto_curated_subset_anndata)
    
    velocyto_curated_subset_anndata.uns['manual_l4_colors']=['#F07E1A', '#B2D06C','#F9CDE2','#F8B366','#BB80B6', "#E31D1E"] # TX monocytes macrophages
    # velocyto_curated_subset_anndata.uns['manual_l4_colors']=['#B2D06C','#F9CDE2','#F8B366','#BB80B6'] # HC macrophages
    # scv.pl.velocity_embedding_stream(velocyto_curated_subset_anndata, basis="umap", color="manual_l4", save = "hc_pf_macrophages.pdf", legend_loc="best", dpi=300, figsize=[6.5,6.5], alpha = 1, density = 1, size = 200)
    scv.pl.velocity_embedding_stream(velocyto_curated_subset_anndata, basis="umap", color="manual_l4", save = "crcpmp_tx_monocytes_macrophages.pdf", legend_loc="best", dpi=300, figsize=[6.5,6.5], alpha = 1, density = 1, size = 200)
    
    velocyto_curated_subset_anndata.write_h5ad(scvelo_anndata_h5ad)
    velocyto_curated_subset_anndata.obs.to_csv(scvelo_anndata_cellmetadata_csv)
    
    return
    
if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--cellmetadata_csv", dest="cellmetadata_csv", required=True)
    argpar.add_argument("--counts_mtx", dest="counts_mtx", required=True)
    argpar.add_argument("--features_csv", dest="features_csv", required=True)
    argpar.add_argument("--umap_csv", dest="umap_csv", required=True)
    argpar.add_argument("--pca_csv", dest="pca_csv", required=True)
    argpar.add_argument("--velocyto_curated_loom", dest="velocyto_curated_loom", required=True)
    argpar.add_argument("--scvelo_anndata_h5ad", dest="scvelo_anndata_h5ad", required=True)
    argpar.add_argument("--scvelo_anndata_cellmetadata_csv", dest="scvelo_anndata_cellmetadata_csv", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.cellmetadata_csv,
        parse_results.counts_mtx,
        parse_results.features_csv,
        parse_results.umap_csv,
        parse_results.pca_csv,
        parse_results.velocyto_curated_loom,
        parse_results.scvelo_anndata_h5ad,
        parse_results.scvelo_anndata_cellmetadata_csv
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)
