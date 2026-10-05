"""
download_gtex.py — CLAUD-IA WP1: GTEx v10 normal-tissue baseline.

Usage:
    python scripts/download_gtex.py

Downloads GTEx v10 bulk RNA-seq and filters to GI-relevant tissues.
These serve as normal-tissue baseline for tumour-vs-normal comparisons,
especially where TCGA has few paired normals (PAAD: 4, CHOL: 9, ESCA: 13).

Output:
    dataset/GTEx/
    ├── raw/
    │   ├── gene_reads.gct.gz   (~898 MB)
    │   ├── sample_attributes.txt
    │   └── subject_phenotypes.txt
    └── processed/
        ├── counts_matrix_GI.csv
        ├── sample_metadata_GI.csv
        └── tissue_summary.txt
"""

import os
import sys
import gzip
import subprocess
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, GTEX_TISSUE_MAP
from run_logger import start_logging

BASE_DIR = os.path.join(DATASET_DIR, "GTEx")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUT_DIR = os.path.join(BASE_DIR, "processed")

os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

URLS = {
    "gene_reads.gct.gz":
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
        "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz",
    "sample_attributes.txt":
        "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/"
        "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
    "subject_phenotypes.txt":
        "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/"
        "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
}


def download_file(fname, url):
    out_path = os.path.join(RAW_DIR, fname)
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        size_mb = os.path.getsize(out_path) / (1024 * 1024)
        print(f"    {fname} exists ({size_mb:.1f} MB) — skip")
        return
    print(f"    Downloading {fname}...")
    result = subprocess.run(["wget", "-q", "--show-progress", "-O", out_path, url])
    if result.returncode != 0:
        print(f"    ERROR downloading {fname}")
        sys.exit(1)


def main():
    print(f"\n{'='*60}")
    print("  CLAUD-IA WP1 — GTEx v10 Normal Tissue Baseline")
    print(f"{'='*60}")

    # 1. Download
    print("\n  [1/4] Downloading GTEx v10 files...")
    for fname, url in URLS.items():
        download_file(fname, url)

    # 2. Load sample annotations
    print("\n  [2/4] Reading sample annotations...")
    attr = pd.read_csv(os.path.join(RAW_DIR, "sample_attributes.txt"),
                       sep="\t", low_memory=False)
    print(f"    Total GTEx samples: {len(attr)}")

    # 3. Filter to GI tissues
    print("\n  [3/4] Filtering to GI-relevant tissues...")
    gi = attr[attr["SMTSD"].isin(GTEX_TISSUE_MAP.keys())].copy()
    gi["WP1_label"] = gi["SMTSD"].map(GTEX_TISSUE_MAP)
    print(f"    GI samples: {len(gi)}")
    for label, count in gi["WP1_label"].value_counts().items():
        print(f"      {count:4d}  {label}")

    gi_ids = set(gi["SAMPID"].tolist())

    # 4. Extract GI columns from counts matrix
    print("\n  [4/4] Extracting GI columns from gene_reads matrix...")
    reads_path = os.path.join(RAW_DIR, "gene_reads.gct.gz")

    with gzip.open(reads_path, "rt") as f:
        f.readline()  # version
        f.readline()  # dimensions
        header = f.readline().strip().split("\t")

    sample_cols = header[2:]
    keep_idx = [0, 1] + [i for i, s in enumerate(sample_cols, start=2) if s in gi_ids]
    print(f"    Matched {len(keep_idx)-2} GI samples in matrix")

    df = pd.read_csv(reads_path, sep="\t", skiprows=2, usecols=keep_idx,
                     compression="gzip", dtype={0: str, 1: str})
    print(f"    Matrix: {df.shape[0]} genes x {df.shape[1]-2} samples")

    # Save
    df.to_csv(os.path.join(OUT_DIR, "counts_matrix_GI.csv"), index=False)
    gi.to_csv(os.path.join(OUT_DIR, "sample_metadata_GI.csv"), index=False)

    with open(os.path.join(OUT_DIR, "tissue_summary.txt"), "w") as f:
        f.write("GTEx v10 GI-tissue sample counts\n" + "="*50 + "\n\n")
        for tissue, label in GTEX_TISSUE_MAP.items():
            n = (gi["SMTSD"] == tissue).sum()
            f.write(f"  {n:4d}  {tissue:50s} -> {label}\n")
        f.write(f"\n  Total: {len(gi)} samples\n")

    # CLDN1 sense check
    cldn1 = df[df.iloc[:, 1] == "CLDN1"]
    if not cldn1.empty:
        vals = cldn1.iloc[0, 2:]
        print(f"    CLDN1: detected in {(vals > 0).sum()}/{len(vals)} samples, "
              f"median={vals.median():.0f}")

    print(f"\n  Done: {BASE_DIR}")


if __name__ == "__main__":
    start_logging()
    main()
