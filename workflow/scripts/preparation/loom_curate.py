#!/usr/bin/env python3
# This script aims to annotate the combined velocyto .loom file.

import os
import shutil
import argparse
import logging
import sys
import re
import loompy
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)

def main(merged_loom, cellmetadata_csv, annotated_loom):
    logging.info("Creating the curated loom file")
    shutil.copy2(merged_loom, annotated_loom)
    logging.info("Importing the annotation file")
    
    cellmetadata = pd.read_csv(cellmetadata_csv)
    
    logging.info("Importing the loom file")
    with loompy.connect(annotated_loom) as ds:
      ds.ca['CellID'] = [re.sub("(.+):([ATCG]{16})x", "S\\1_\\2", x) for x in ds.ca['CellID']]
      ds.ca['CellID'] = [re.sub("-", "_", x) for x in ds.ca['CellID']]
      ds.ca['CellID'] = [x + "-1" for x in ds.ca['CellID']]
      
      ca_pd = pd.DataFrame(ds.ca['CellID'], columns = ['CellID'])
      ca_pd = pd.merge(ca_pd, cellmetadata, how='left', on = "CellID")
      
      ds.ca['SampleID'] = ca_pd['SampleID'].to_numpy()
      ds.ca['Donor'] = ca_pd['Donor'].to_numpy()
      ds.ca['Group'] = ca_pd['Group'].to_numpy()
      ds.ca['Tissue'] = ca_pd['Tissue'].to_numpy()
      ds.ca['manual_l1'] = ca_pd['manual_l1'].to_numpy()
      ds.ca['manual_l2'] = ca_pd['manual_l2'].to_numpy()
      ds.ca['manual_l3'] = ca_pd['manual_l3'].to_numpy()
      ds.ca['manual_l4'] = ca_pd['manual_l4'].to_numpy()
    
    logging.info("Finished combining loom file")
    return
    
if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--merged_loom", dest="merged_loom", required=True)
    argpar.add_argument("--cellmetadata_csv", dest="cellmetadata_csv", required=True)
    argpar.add_argument("--annotated_loom", dest="annotated_loom", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.merged_loom,
        parse_results.cellmetadata_csv,
        parse_results.annotated_loom,
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)
