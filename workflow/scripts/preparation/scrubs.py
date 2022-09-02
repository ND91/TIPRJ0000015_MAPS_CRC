import scrublet
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse
import logging
import sys
import re

logging.basicConfig(
    level=logging.INFO, format="%(name)s - %(levelname)s %(asctime)s - %(message)s"
)


def main(matrix_path, genes_path, barcode_path, output, sampleid):
    logging.info("Generating paths")
    # create output directory and ignore error if dir is already there
    try:
        os.makedirs(output)
    except FileExistsError:
        if os.path.isdir(output):
            pass
        else:
            os.error("output file already exists and is not a directory")

    logging.info("reading matrix")
    counts_matrix = scipy.io.mmread(matrix_path).T.tocsc()
    logging.info("loading genes")
    genes = np.array(scrublet.load_genes(genes_path, delimiter="\t", column=1))
    logging.info("loading barcode")
    cb_dirty = np.loadtxt(barcode_path, dtype="str")
    cb = np.array([re.sub("-[0-9]+", "_" + sampleid, x) for x in cb_dirty])
    logging.info("Creating scrublet object")
    scrub = scrublet.Scrublet(counts_matrix, expected_doublet_rate=0.08)
    logging.info("Calculating doublets")
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
    )
    logging.info("creating histogram")
    hist_fig, hist_axs = scrub.plot_histogram()
    logging.info("saving histogram")
    hist_fig.savefig(f"{output}/{sampleid}_hist_scrublet.png")
    logging.info("set embedded cells")
    scrub.set_embedding(
        "UMAP", scrublet.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
    )
    logging.info("creating UMAP plot")
    fig, axs = scrub.plot_embedding("UMAP", order_points=True)
    logging.info("save UMAP plot")
    fig.savefig(f"{output}/{sampleid}_scatterplot_scrublet.png")
    logging.info("creating scrublet summary table")
    scrublet_summary = pd.DataFrame(
        {
            "CB": cb,
            "Doublet scores": doublet_scores,
            "Predicted doublets": predicted_doublets,
        }
    )
    logging.info("save scrublet summary table to csv.gz file")
    scrublet_summary.to_csv(
        f"{output}/{sampleid}_scrublet.csv.gz", index=False, compression="gzip"
    )
    logging.info(" done save scrublet summary table to csv.gz file")
    logging.info("all done")
    return


if __name__ == "__main__":
    argpar = argparse.ArgumentParser()
    argpar.add_argument("--matrix", dest="matrix", required=True)
    argpar.add_argument("--genes", dest="genes", required=True)
    argpar.add_argument("--barcode", dest="barcode", required=True)
    argpar.add_argument("--output", dest="output", required=True)
    argpar.add_argument("--sampleid", dest="sampleid", required=True)

    parse_results = argpar.parse_args()

    main(
        parse_results.matrix,
        parse_results.genes,
        parse_results.barcode,
        parse_results.output,
        parse_results.sampleid,
    )
    # exiting this hard way is need because a gentle exit keeps the script hanging on workernodes
    os._exit(0)
