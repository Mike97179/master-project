"""
summary_report.py — CLAUD-IA WP1: dataset summary table for reporting.

Scans whatever has been downloaded under dataset/TCGA/ and writes one row
per cancer type to dataset/summary/tcga_summary.{csv,xlsx}: sample counts
by type, excluded samples, matrix dimensions and CLDN1 statistics.

Usage:
    python scripts/summary_report.py
    python scripts/summary_report.py --md5     # also verify checksums (slow)
"""

import os
import sys
import json
import hashlib

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES, TCGA_VALID_SAMPLE_CODES
from run_logger import start_logging

SUMMARY_DIR = os.path.join(DATASET_DIR, "summary")


def cldn1_stats(matrix_path):
    """Pull the CLDN1 row without loading the whole matrix into memory.

    The matrices run to ~80 MB each; a plain read_csv of all six would cost
    far more than a single line scan per file.
    """
    if not os.path.exists(matrix_path):
        return {}
    with open(matrix_path) as f:
        f.readline()                      # header: gene_id,gene_name,gene_type,<barcodes>
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) > 2 and parts[1] == "CLDN1":
                vals = pd.to_numeric(pd.Series(parts[3:]), errors="coerce").dropna()
                return {
                    "cldn1_detected": int((vals > 0).sum()),
                    "cldn1_median": round(float(vals.median()), 1),
                    "cldn1_min": int(vals.min()),
                    "cldn1_max": int(vals.max()),
                }
    return {}


def verify_md5(rna_dir):
    """Count files whose checksum matches the manifest."""
    manifest = os.path.join(rna_dir, "manifest.txt")
    if not os.path.exists(manifest):
        return None
    ok = 0
    rows = [l.rstrip("\n").split("\t") for l in open(manifest)][1:]
    for fid, fname, md5, *_ in rows:
        path = os.path.join(rna_dir, "raw", fid, fname)
        if not os.path.isfile(path):
            continue
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() == md5:
            ok += 1
    return ok


def dir_size_gb(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return round(total / (1024 ** 3), 2)


def collect(project_id, check_md5=False):
    base = os.path.join(DATASET_DIR, "TCGA", project_id)
    rna_dir = os.path.join(base, "rna")
    processed = os.path.join(base, "processed")

    row = {
        "project": project_id,
        "cancer_type": CANCER_TYPES[project_id]["label"],
        "status": "not downloaded",
    }

    if not os.path.isdir(base):
        return row

    # Files available at the GDC, as recorded at download time
    meta_path = os.path.join(rna_dir, "metadata.json")
    if os.path.exists(meta_path):
        row["files_at_gdc"] = len(json.load(open(meta_path)))

    raw_dir = os.path.join(rna_dir, "raw")
    if os.path.isdir(raw_dir):
        row["files_downloaded"] = sum(
            1 for r, _, fs in os.walk(raw_dir)
            for f in fs if f.endswith(".tsv") and not f.endswith(".partial")
        )
        row["files_incomplete"] = sum(
            1 for r, _, fs in os.walk(raw_dir) for f in fs if f.endswith(".partial")
        )

    if check_md5:
        row["md5_verified"] = verify_md5(rna_dir)

    # Sample breakdown
    meta_csv = os.path.join(processed, "sample_metadata.csv")
    if os.path.exists(meta_csv):
        sm = pd.read_csv(meta_csv)
        row["samples_total"] = len(sm)
        row["patients_unique"] = sm["patient"].nunique()
        counts = sm["sample_type_label"].value_counts()
        row["primary_tumor"] = int(counts.get("Primary Tumor", 0))
        row["solid_tissue_normal"] = int(counts.get("Solid Tissue Normal", 0))
        row["metastatic"] = int(counts.get("Metastatic", 0))

    # Excluded samples, with the codes spelled out for the write-up
    excl_csv = os.path.join(processed, "excluded_samples.csv")
    if os.path.exists(excl_csv):
        ex = pd.read_csv(excl_csv, dtype={"sample_type_code": str})
        row["excluded_total"] = len(ex)
        if len(ex):
            row["excluded_detail"] = "; ".join(
                f"{c} {l} x{n}" for (c, l), n in
                ex.groupby(["sample_type_code", "sample_type_label"]).size().items()
            )

    matrix_path = os.path.join(processed, "counts_matrix.csv")
    if os.path.exists(matrix_path):
        with open(matrix_path) as f:
            row["genes"] = sum(1 for _ in f) - 1
        row.update(cldn1_stats(matrix_path))

    row["clinical"] = "yes" if os.path.exists(
        os.path.join(base, "clinical", "clinical.tsv")) else "no"
    row["size_gb"] = dir_size_gb(base)

    expected = row.get("files_at_gdc")
    got = row.get("files_downloaded", 0)
    if expected and got >= expected and os.path.exists(matrix_path):
        row["status"] = "complete"
    elif got:
        row["status"] = f"partial ({got}/{expected})"
    return row


COLUMNS = [
    "project", "cancer_type", "status",
    "files_at_gdc", "files_downloaded", "files_incomplete", "md5_verified",
    "samples_total", "primary_tumor", "solid_tissue_normal", "metastatic",
    "patients_unique", "excluded_total", "excluded_detail",
    "genes", "cldn1_detected", "cldn1_median", "cldn1_min", "cldn1_max",
    "clinical", "size_gb",
]


def main():
    check_md5 = "--md5" in sys.argv
    print(f"\n  Sample codes kept: {', '.join(sorted(TCGA_VALID_SAMPLE_CODES))}")
    if check_md5:
        print("  Verifying md5 checksums (this takes a few minutes)...")

    rows = []
    for project_id in CANCER_TYPES:
        print(f"  Reading {project_id}...")
        rows.append(collect(project_id, check_md5))

    df = pd.DataFrame(rows)
    df = df.reindex(columns=[c for c in COLUMNS if c in df.columns])

    # Totals row, over the numeric columns only
    totals = {"project": "TOTAL", "cancer_type": "", "status": ""}
    for col in df.columns:
        if col in totals or not pd.api.types.is_numeric_dtype(df[col]):
            continue
        totals[col] = df[col].sum()
    # Columns where a sum is meaningless: the same gene list repeated six
    # times, and statistics that do not add up.
    for col in ["genes", "cldn1_median", "cldn1_min", "cldn1_max"]:
        totals.pop(col, None)
    df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)

    # Nullable ints: a column with gaps would otherwise render as 523.0.
    # Columns with real decimals (size_gb, cldn1_median) are left alone.
    for col in df.columns:
        if not pd.api.types.is_numeric_dtype(df[col]):
            continue
        vals = df[col].dropna()
        if len(vals) and (vals % 1 == 0).all():
            df[col] = df[col].astype("Int64")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    csv_path = os.path.join(SUMMARY_DIR, "tcga_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n  CSV:   {os.path.relpath(csv_path, os.path.dirname(DATASET_DIR))}")

    try:
        xlsx_path = os.path.join(SUMMARY_DIR, "tcga_summary.xlsx")
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="TCGA summary", index=False)
            sheet = writer.sheets["TCGA summary"]
            for i, col in enumerate(df.columns, start=1):
                cells = ["" if pd.isna(v) else str(v) for v in df[col]]
                width = max(len(col), *(len(c) for c in cells)) + 2
                sheet.column_dimensions[sheet.cell(1, i).column_letter].width = min(width, 50)
            sheet.freeze_panes = "A2"
        print(f"  Excel: {os.path.relpath(xlsx_path, os.path.dirname(DATASET_DIR))}")
    except ImportError:
        print("  Excel skipped (pip install openpyxl)")

    print()
    cols = ["project", "status", "samples_total", "primary_tumor",
            "solid_tissue_normal", "metastatic", "excluded_total", "cldn1_median"]
    view = df[[c for c in cols if c in df.columns]]
    view = view.apply(lambda c: ["-" if pd.isna(v) else str(v) for v in c])
    print(view.to_string(index=False))


if __name__ == "__main__":
    start_logging()
    main()
