"""
download_gtex.py — CLAUD-IA WP1: GTEx v10 normal-tissue baseline.

Usage:
    python scripts/download_gtex.py
    python scripts/download_gtex.py --verify     # check files, download nothing

Downloads GTEx v10 bulk RNA-seq and filters to GI-relevant tissues.
These serve as normal-tissue baseline for tumour-vs-normal comparisons,
especially where TCGA has few paired normals (PAAD: 4, CHOL: 9, ESCA: 13).

Sample filter:
    - GI tissues only (see GTEX_TISSUE_MAP in config.py)
    - SMAFRZE == RNASEQ: of the 8223 GI samples, only 3714 were bulk
      RNA-sequenced. The rest are miRNA (SMLRNA), WGS/WES, or samples
      GTEx itself marked EXCLUDE, and none of them appear in the counts
      matrix. Keeping them would leave the metadata out of step with the
      matrix columns.

Output:
    dataset/GTEx/
    ├── raw/
    │   ├── gene_reads.gct.gz   (~898 MB)
    │   ├── sample_attributes.txt
    │   └── subject_phenotypes.txt
    └── processed/
        ├── counts_matrix_GI.csv
        ├── sample_metadata_GI.csv
        ├── excluded_samples.csv
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

# Rough expected sizes, to catch a truncated download before it is parsed
MIN_SIZES = {
    "gene_reads.gct.gz": 800 * 1024 * 1024,
    "sample_attributes.txt": 30 * 1024 * 1024,
    "subject_phenotypes.txt": 10 * 1024,
}


def download_file(fname, url):
    out_path = os.path.join(RAW_DIR, fname)
    if os.path.exists(out_path) and os.path.getsize(out_path) >= MIN_SIZES[fname]:
        size_mb = os.path.getsize(out_path) / (1024 * 1024)
        print(f"    {fname} exists ({size_mb:.1f} MB) — skip")
        return
    print(f"    Downloading {fname}...")
    # -c resumes a partial file instead of starting over
    result = subprocess.run(["wget", "-c", "-q", "--show-progress", "-O", out_path, url])
    if result.returncode != 0:
        print(f"    ERROR downloading {fname}")
        sys.exit(1)


def verify_files():
    """Check every raw file is present, full-sized and not corrupt.

    GTEx publishes no checksums, unlike the GDC. A gzip file carries its
    own CRC32, so decompressing it end to end is the equivalent integrity
    check; the plain text files are only size-checked.
    """
    ok = True
    for fname in URLS:
        path = os.path.join(RAW_DIR, fname)
        if not os.path.exists(path):
            print(f"    {fname}: MISSING")
            ok = False
            continue
        size = os.path.getsize(path)
        if size < MIN_SIZES[fname]:
            print(f"    {fname}: TOO SMALL ({size/1048576:.1f} MB) — truncated?")
            ok = False
            continue
        if fname.endswith(".gz"):
            print(f"    {fname}: verifying gzip CRC ({size/1048576:.1f} MB)...")
            if subprocess.run(["gzip", "-t", path]).returncode != 0:
                print(f"    {fname}: CORRUPT (CRC mismatch)")
                ok = False
                continue
        print(f"    {fname}: OK ({size/1048576:.1f} MB)")
    return ok


def load_gi_samples():
    """GI samples that were actually bulk RNA-sequenced, plus the rest."""
    attr = pd.read_csv(os.path.join(RAW_DIR, "sample_attributes.txt"),
                       sep="\t", low_memory=False)
    print(f"    Total GTEx samples (all assays): {len(attr)}")

    gi = attr[attr["SMTSD"].isin(GTEX_TISSUE_MAP)].copy()
    gi["WP1_label"] = gi["SMTSD"].map(GTEX_TISSUE_MAP)
    print(f"    In GI tissues: {len(gi)}")

    rna = gi[gi["SMAFRZE"] == "RNASEQ"].copy()
    other = gi[gi["SMAFRZE"] != "RNASEQ"].copy()
    print(f"    Bulk RNA-seq (SMAFRZE=RNASEQ): {len(rna)}")
    if len(other):
        print(f"    Excluded {len(other)} GI samples from other assays:")
        for assay, n in other["SMAFRZE"].value_counts(dropna=False).items():
            print(f"      {str(assay):10} {n}")

    # Donor sex / age group / death circumstances
    pheno = pd.read_csv(os.path.join(RAW_DIR, "subject_phenotypes.txt"), sep="\t")
    rna["SUBJID"] = rna["SAMPID"].str.split("-").str[:2].str.join("-")
    rna = rna.merge(pheno, on="SUBJID", how="left")
    print(f"    Donor phenotypes merged: {rna['SEX'].notna().sum()}/{len(rna)} rows")
    return rna, other


def extract_matrix(rna):
    """Pull just the GI sample columns out of the 898 MB counts matrix."""
    reads_path = os.path.join(RAW_DIR, "gene_reads.gct.gz")
    with gzip.open(reads_path, "rt") as f:
        f.readline()                                  # '#1.2'
        n_genes, n_samples = f.readline().split("\t")  # dimensions
        header = f.readline().rstrip("\n").split("\t")
    print(f"    Matrix on disk: {n_genes.strip()} genes x {n_samples.strip()} samples")

    wanted = set(rna["SAMPID"])
    keep = [0, 1] + [i for i, s in enumerate(header) if i >= 2 and s in wanted]
    print(f"    Matched {len(keep) - 2}/{len(wanted)} GI samples as columns")

    # int32 keeps the frame near 900 MB instead of 1.8 GB; GTEx read
    # counts stay far below the int32 ceiling.
    names = [header[i] for i in keep]
    dtypes = {n: "int32" for n in names[2:]}
    dtypes[names[0]] = str
    dtypes[names[1]] = str

    df = pd.read_csv(reads_path, sep="\t", skiprows=2, usecols=keep, dtype=dtypes)
    print(f"    Extracted: {df.shape[0]} genes x {df.shape[1] - 2} samples")
    return df


def write_outputs(df, rna, other):
    os.makedirs(OUT_DIR, exist_ok=True)

    # Keep the metadata in step with the matrix columns, in the same order.
    # GTEx happens to ship sample_attributes.txt in the .gct column order, but
    # that is luck, not a guarantee: reindexing makes the correspondence
    # explicit so a positional pairing cannot silently go wrong.
    col_order = list(df.columns[2:])
    rna = (rna[rna["SAMPID"].isin(set(col_order))]
           .set_index("SAMPID")
           .reindex(col_order)
           .reset_index())
    assert list(rna["SAMPID"]) == col_order

    df.to_csv(os.path.join(OUT_DIR, "counts_matrix_GI.csv"), index=False)
    rna.to_csv(os.path.join(OUT_DIR, "sample_metadata_GI.csv"), index=False)

    cols = [c for c in ["SAMPID", "SMTSD", "WP1_label", "SMAFRZE"] if c in other.columns]
    other[cols].to_csv(os.path.join(OUT_DIR, "excluded_samples.csv"), index=False)

    with open(os.path.join(OUT_DIR, "tissue_summary.txt"), "w") as f:
        f.write("GTEx v10 GI-tissue bulk RNA-seq sample counts\n" + "=" * 60 + "\n\n")
        for tissue, label in GTEX_TISSUE_MAP.items():
            n = (rna["SMTSD"] == tissue).sum()
            f.write(f"  {n:5d}  {tissue:48s} -> {label}\n")
        f.write(f"\n  Total: {len(rna)} samples\n\n")
        f.write("Grouped by WP1 cancer type:\n")
        for label, n in rna["WP1_label"].value_counts().items():
            f.write(f"  {n:5d}  {label}\n")

    print(f"    counts_matrix_GI.csv:   {df.shape[0]} x {df.shape[1] - 2}")
    print(f"    sample_metadata_GI.csv: {len(rna)} rows, aligned to matrix columns")
    print(f"    excluded_samples.csv:   {len(other)} rows")
    return rna


def sense_check(df, rna):
    cldn1 = df[df["Description"] == "CLDN1"]
    if cldn1.empty:
        print("    WARNING: CLDN1 not found")
        return
    vals = cldn1.iloc[0, 2:].astype("int64")
    print(f"    CLDN1: detected in {(vals > 0).sum()}/{len(vals)} samples, "
          f"median={vals.median():.0f}, range=[{vals.min()}, {vals.max()}]")
    print("\n    CLDN1 median by WP1 group (raw counts, not yet normalised):")
    label_of = dict(zip(rna["SAMPID"], rna["WP1_label"]))
    by_group = vals.groupby(vals.index.map(label_of)).median().sort_values()
    for label, med in by_group.items():
        print(f"      {label:12} {med:10.0f}")


def main():
    verify_only = "--verify" in sys.argv
    print(f"\n{'=' * 60}")
    print("  CLAUD-IA WP1 — GTEx v10 Normal Tissue Baseline")
    print(f"{'=' * 60}")

    os.makedirs(RAW_DIR, exist_ok=True)

    if verify_only:
        print("\n  Verifying raw files (no download)...")
        sys.exit(0 if verify_files() else 1)

    print("\n  [1/5] Downloading GTEx v10 files...")
    for fname, url in URLS.items():
        download_file(fname, url)

    print("\n  [2/5] Verifying downloaded files...")
    if not verify_files():
        print("    Integrity check failed — re-run to resume the download")
        sys.exit(1)

    print("\n  [3/5] Selecting GI bulk RNA-seq samples...")
    rna, other = load_gi_samples()

    print("\n  [4/5] Extracting GI columns from the counts matrix...")
    df = extract_matrix(rna)
    rna = write_outputs(df, rna, other)

    print("\n  [5/5] Sense check...")
    sense_check(df, rna)

    print(f"\n  Done: {BASE_DIR}")


if __name__ == "__main__":
    start_logging()
    main()
