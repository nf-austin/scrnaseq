#!/usr/bin/env python3

import argparse
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.stats import median_abs_deviation
import warnings

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

# Species-specific gene prefix conventions.
# Human Ensembl/HGNC: uppercase (MT-, RPS/RPL, HB*).
# Mouse Ensembl/MGI: capitalized (mt-, Rps/Rpl, Hb*).
SPECIES_PREFIXES = {
    "human": {"mt": ("MT-",),        "ribo": ("RPS", "RPL"), "hb": r"^HB[^P]"},
    "mouse": {"mt": ("mt-",),        "ribo": ("Rps", "Rpl"), "hb": r"^Hb[^p]"},
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5', required=True)
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--species', required=True, choices=sorted(SPECIES_PREFIXES.keys()))
    parser.add_argument('--run_scrublet', action='store_true')
    args = parser.parse_args()

    prefixes = SPECIES_PREFIXES[args.species]

    # Load Data
    adata = sc.read_10x_h5(args.h5)
    adata.var_names_make_unique()

    # Annotate gene categories (species-aware)
    adata.var["mt"] = adata.var_names.str.startswith(prefixes["mt"])
    adata.var["ribo"] = adata.var_names.str.startswith(prefixes["ribo"])
    adata.var["hb"] = adata.var_names.str.contains(prefixes["hb"])

    # Calculate metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        percent_top=None,
        log1p=True
    )

    # MAD-based filtering logic
    def is_outlier(adata, metric, nmads):
        M = adata.obs[metric]
        mad = median_abs_deviation(M)
        if mad == 0:
            return np.zeros(len(M), dtype=bool)
        outlier = (M < np.median(M) - nmads * mad) | (M > np.median(M) + nmads * mad)
        return outlier

    adata.obs["outlier_counts"] = is_outlier(adata, "log1p_total_counts", 5)
    adata.obs["outlier_genes"] = is_outlier(adata, "log1p_n_genes_by_counts", 5)
    adata.obs["outlier_mt"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8.0)

    # Base mask
    base_mask = (
        adata.obs["outlier_counts"] |
        adata.obs["outlier_genes"] |
        adata.obs["outlier_mt"]
    )

    if args.run_scrublet:
        sc.pp.scrublet(adata)
        adata.obs["passing_qc"] = ~(base_mask | adata.obs.get("predicted_doublet", False))
    else:
        adata.obs["passing_qc"] = ~base_mask

    # Save to h5ad
    adata.write_h5ad(f"{args.sample_id}_annotated.h5ad")

if __name__ == "__main__":
    main()