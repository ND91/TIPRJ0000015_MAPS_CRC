#!/usr/bin/env python3
# This script aims to combine the velocyto .loom files into a single loom file.

import os
import argparse
import logging
import sys
import re
import loompy

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)

def main(velocyto_looms, combined_loom):
    logging.info("Combining loom files")
    
    loompy.combine(files = velocyto_looms, 
                   key = "Accession", 
                   output_file = f"{combined_loom}")
    logging.info("Finished combining loom file")
    return
    
if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--velocyto_looms", dest="velocyto_looms", required=True, nargs='+')
    argpar.add_argument("--combined_loom", dest="combined_loom", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.velocyto_looms,
        parse_results.combined_loom,
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)
