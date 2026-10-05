"""
inventory.py — CLAUD-IA WP1: Show what's been collected so far.

Usage:
    python scripts/inventory.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES
from run_logger import start_logging

def count_files(path, ext=None):
    if not os.path.isdir(path):
        return 0
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            if ext is None or f.endswith(ext):
                total += 1
    return total

def file_size_mb(path):
    if os.path.isfile(path):
        return os.path.getsize(path) / (1024 * 1024)
    return 0

def main():
    print(f"\n{'='*70}")
    print("  CLAUD-IA WP1 — Dataset Inventory")
    print(f"{'='*70}\n")

    # TCGA
    print("  TCGA (bulk RNA-seq)")
    print("  " + "-"*66)
    tcga_dir = os.path.join(DATASET_DIR, "TCGA")
    for project_id, info in CANCER_TYPES.items():
        base = os.path.join(tcga_dir, project_id)
        raw_dir = os.path.join(base, "rna", "raw")
        counts = os.path.join(base, "processed", "counts_matrix.csv")
        clinical = os.path.join(base, "clinical", "clinical.tsv")

        raw_n = count_files(raw_dir, ".tsv")
        has_counts = "yes" if os.path.isfile(counts) else "no"
        has_clin = "yes" if os.path.isfile(clinical) else "no"

        print(f"    {project_id:12s}  raw_files={raw_n:4d}  "
              f"counts={'yes' if os.path.isfile(counts) else 'no':3s}  "
              f"clinical={'yes' if os.path.isfile(clinical) else 'no':3s}")

    # GTEx
    print(f"\n  GTEx (normal tissue baseline)")
    print("  " + "-"*66)
    gtex_counts = os.path.join(DATASET_DIR, "GTEx", "processed", "counts_matrix_GI.csv")
    if os.path.isfile(gtex_counts):
        print(f"    counts_matrix_GI.csv: {file_size_mb(gtex_counts):.1f} MB")
    else:
        print("    Not downloaded yet")

    # GEO
    print(f"\n  GEO (bulk RNA-seq datasets)")
    print("  " + "-"*66)
    geo_dir = os.path.join(DATASET_DIR, "GEO")
    if os.path.isdir(geo_dir):
        for entry in sorted(os.listdir(geo_dir)):
            entry_path = os.path.join(geo_dir, entry)
            if os.path.isdir(entry_path):
                expr = os.path.join(entry_path, "expression.tsv")
                has_expr = f"{file_size_mb(expr):.1f} MB" if os.path.isfile(expr) else "no"
                print(f"    {entry:30s}  expression={has_expr}")
        search = os.path.join(geo_dir, "search_results.csv")
        if os.path.isfile(search):
            print(f"    search_results.csv present")
    else:
        print("    No GEO data yet")

    # PubMed search results
    pubmed_csv = os.path.join(DATASET_DIR, "pubmed_datasets.csv")
    if os.path.isfile(pubmed_csv):
        import pandas as pd
        df = pd.read_csv(pubmed_csv)
        print(f"\n  PubMed dataset search: {len(df)} accessions found")

    # GSA/NODE search results
    gsa_csv = os.path.join(DATASET_DIR, "gsa_node_results.csv")
    if os.path.isfile(gsa_csv):
        import pandas as pd
        df = pd.read_csv(gsa_csv)
        print(f"\n  GSA/NODE search: {len(df)} datasets found")

    # Disk usage
    print(f"\n  Total disk usage")
    print("  " + "-"*66)
    total_size = 0
    for root, _, files in os.walk(DATASET_DIR):
        for f in files:
            total_size += os.path.getsize(os.path.join(root, f))
    print(f"    {DATASET_DIR}: {total_size / (1024**3):.2f} GB")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    start_logging()
    main()
