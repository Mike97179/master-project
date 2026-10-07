"""
summary_report.py — CLAUD-IA WP1: dataset summary table for reporting.

Scans whatever has been downloaded under dataset/ and writes, into
dataset/summary/:

    tcga_summary.csv   one row per TCGA cancer type
    gtex_summary.csv   one row per GTEx GI tissue
    wp1_summary.xlsx   both of the above plus "Normals per cancer", which
                       is the table the project hinges on: how many normal
                       samples each cancer type ends up with once GTEx is
                       added to TCGA's few paired normals.

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
from config import (DATASET_DIR, CANCER_TYPES, TCGA_VALID_SAMPLE_CODES,
                    GTEX_TISSUE_MAP)
from run_logger import start_logging

SUMMARY_DIR = os.path.join(DATASET_DIR, "summary")


def check_alignment(matrix_path, meta_path, id_col, n_lead):
    """Do the metadata rows line up with the matrix columns, in order?

    Same length and same ids is not enough: a positional pairing needs the
    same order too, and the two tables are built from unrelated orderings.
    """
    if not (os.path.exists(matrix_path) and os.path.exists(meta_path)):
        return None
    with open(matrix_path) as f:
        cols = f.readline().rstrip("\n").split(",")[n_lead:]
    rows = list(pd.read_csv(meta_path, usecols=[id_col], low_memory=False)[id_col])
    if set(cols) != set(rows):
        return "id mismatch"
    return "yes" if cols == rows else "order differs"


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

    row["meta_aligned"] = check_alignment(matrix_path, meta_csv, "barcode", 3)
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


def cldn1_stats_gtex(matrix_path, sampids=None):
    """CLDN1 row of the GTEx matrix; 'Description' holds the gene symbol."""
    if not os.path.exists(matrix_path):
        return {}
    with open(matrix_path) as f:
        header = f.readline().rstrip("\n").split(",")
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) > 2 and parts[1] == "CLDN1":
                vals = pd.Series(
                    pd.to_numeric(pd.Series(parts[2:]), errors="coerce").values,
                    index=header[2:],
                ).dropna()
                if sampids is not None:
                    vals = vals[vals.index.isin(sampids)]
                if vals.empty:
                    return {}
                return {
                    "cldn1_detected": int((vals > 0).sum()),
                    "cldn1_median": round(float(vals.median()), 1),
                    "cldn1_min": int(vals.min()),
                    "cldn1_max": int(vals.max()),
                }
    return {}


def collect_gtex():
    """One row per GI tissue, from the processed GTEx metadata."""
    base = os.path.join(DATASET_DIR, "GTEx")
    meta_csv = os.path.join(base, "processed", "sample_metadata_GI.csv")
    matrix = os.path.join(base, "processed", "counts_matrix_GI.csv")
    if not os.path.exists(meta_csv):
        return pd.DataFrame()

    sm = pd.read_csv(meta_csv, low_memory=False)
    excl_csv = os.path.join(base, "processed", "excluded_samples.csv")
    excl = pd.read_csv(excl_csv) if os.path.exists(excl_csv) else pd.DataFrame()

    rows = []
    for tissue, label in GTEX_TISSUE_MAP.items():
        sub = sm[sm["SMTSD"] == tissue]
        if sub.empty:
            continue
        row = {
            "tissue": tissue,
            "wp1_label": label,
            "samples": len(sub),
            "donors": sub["SUBJID"].nunique() if "SUBJID" in sub else None,
        }
        if "SEX" in sub:
            row["male"] = int((sub["SEX"] == 1).sum())
            row["female"] = int((sub["SEX"] == 2).sum())
        if len(excl) and "SMTSD" in excl:
            row["excluded_other_assays"] = int((excl["SMTSD"] == tissue).sum())
        row.update(cldn1_stats_gtex(matrix, set(sub["SAMPID"])))
        row["meta_aligned"] = check_alignment(matrix, meta_csv, "SAMPID", 2)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    totals = {"tissue": "TOTAL", "wp1_label": ""}
    for col in df.columns:
        if col in totals or not pd.api.types.is_numeric_dtype(df[col]):
            continue
        totals[col] = df[col].sum()
    for col in ["cldn1_median", "cldn1_min", "cldn1_max"]:
        totals.pop(col, None)
    return pd.concat([df, pd.DataFrame([totals])], ignore_index=True)


# GTEx groups liver donors under one label for two cancers, and terminal
# ileum maps to no WP1 cancer at all.
GTEX_LABEL_TO_PROJECTS = {
    "PAAD": ["TCGA-PAAD"],
    "COAD": ["TCGA-COAD"],
    "STAD": ["TCGA-STAD"],
    "ESCA": ["TCGA-ESCA"],
    "LIHC_CHOL": ["TCGA-LIHC", "TCGA-CHOL"],
}


def normals_table(tcga, gtex):
    """Normal samples per cancer type, TCGA vs GTEx.

    This is the point of downloading GTEx: TCGA leaves PAAD with 4 normals
    and CHOL with 9, too few for a tumour-vs-normal contrast.
    """
    if gtex.empty:
        return pd.DataFrame()
    gtex_body = gtex[gtex["tissue"] != "TOTAL"]
    per_project = {}
    for label, projects in GTEX_LABEL_TO_PROJECTS.items():
        n = int(gtex_body.loc[gtex_body["wp1_label"] == label, "samples"].sum())
        for proj in projects:
            per_project[proj] = n

    rows = []
    body = tcga[tcga["project"] != "TOTAL"]
    for _, r in body.iterrows():
        tcga_n = r.get("solid_tissue_normal")
        tcga_n = 0 if pd.isna(tcga_n) else int(tcga_n)
        gtex_n = per_project.get(r["project"], 0)
        rows.append({
            "project": r["project"],
            "cancer_type": r["cancer_type"],
            "tcga_tumour": 0 if pd.isna(r.get("primary_tumor")) else int(r["primary_tumor"]),
            "tcga_normal": tcga_n,
            "gtex_normal": gtex_n,
            "normal_total": tcga_n + gtex_n,
            "note": "GTEx liver shared with CHOL" if r["project"] in ("TCGA-LIHC", "TCGA-CHOL") else "",
        })
    df = pd.DataFrame(rows)
    totals = {"project": "TOTAL", "cancer_type": "", "note": ""}
    for col in ["tcga_tumour", "tcga_normal", "gtex_normal", "normal_total"]:
        totals[col] = df[col].sum()
    return pd.concat([df, pd.DataFrame([totals])], ignore_index=True)


def collect_geo():
    """GEO candidates: what the search found and what inspection concluded.

    GEO is a curation step, not a bulk download: the useful summary is how
    many candidates survive each filter, not sample counts.
    """
    base = os.path.join(DATASET_DIR, "GEO")
    insp_csv = os.path.join(base, "inspected.csv")
    if not os.path.isdir(base):
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    insp = pd.read_csv(insp_csv) if os.path.exists(insp_csv) else pd.DataFrame()

    # Search results and decisions live under dataset/GEO/<CANCER>/
    search_parts, dec_parts, downloaded = [], [], set()
    for short in [v["short"] for v in CANCER_TYPES.values()]:
        cdir = os.path.join(base, short)
        if not os.path.isdir(cdir):
            continue
        for name, bucket in (("search_results.csv", search_parts),
                             ("decisions.csv", dec_parts)):
            path = os.path.join(cdir, name)
            if os.path.exists(path):
                part = pd.read_csv(path)
                part["search_cancer"] = short
                bucket.append(part)
        downloaded |= {d for d in os.listdir(cdir)
                       if d.startswith("GSE") and os.path.isdir(os.path.join(cdir, d))}

    search = pd.concat(search_parts, ignore_index=True) if search_parts else pd.DataFrame()
    decisions = pd.concat(dec_parts, ignore_index=True) if dec_parts else pd.DataFrame()
    if search.empty and insp.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # The funnel the GEO pipeline reports: candidates found, how many
    # survived each gate, and how the survivors were finally classed.
    rows = []
    for short in [v["short"] for v in CANCER_TYPES.values()]:
        sub = search[search["search_cancer"] == short] if len(search) else pd.DataFrame()
        if sub.empty:
            continue
        accs = set(sub["accession"])
        sub_insp = insp[insp["accession"].isin(accs)] if len(insp) else pd.DataFrame()
        sub_dec = (decisions[decisions["accession"].isin(accs)]
                   if len(decisions) else pd.DataFrame())

        row = {"cancer": short, "candidates_found": len(sub)}
        if len(sub_insp):
            hard = sub_insp["hard"].fillna("")
            row["past_gates_1_2"] = int((hard == "").sum())
            row["with_raw_counts"] = int((sub_insp["has_raw"] == "yes").sum())
        if len(sub_dec):
            row["tissue_confirmed"] = int((sub_dec["decision"] == "DOWNLOAD").sum())
            row["ambiguous"] = int((sub_dec["decision"] == "AMBIGUOUS").sum())
            row["rejected"] = int((sub_dec["decision"] == "SKIP").sum())
            if "declares" in sub_dec:
                decl = sub_dec["declares"].fillna("")
                row["declare_tissue"] = int((decl == "tissue").sum())
                row["declare_culture"] = int((decl == "culture").sum())
        row["downloaded"] = len(accs & downloaded)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df, insp, decisions
    totals = {"cancer": "TOTAL"}
    for col in df.columns:
        if col == "cancer":
            continue
        totals[col] = df[col].sum()
    df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
    return df, insp, decisions


def write_excel(path, sheets):
    """One workbook, one sheet per table, columns sized to their contents."""
    try:
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for name, df in sheets.items():
                if df is None or df.empty:
                    continue
                df.to_excel(writer, sheet_name=name, index=False)
                sheet = writer.sheets[name]
                for i, col in enumerate(df.columns, start=1):
                    cells = ["" if pd.isna(v) else str(v) for v in df[col]]
                    width = max(len(col), *(len(c) for c in cells)) + 2
                    sheet.column_dimensions[sheet.cell(1, i).column_letter].width = min(width, 50)
                sheet.freeze_panes = "A2"
        return True
    except ImportError:
        return False


COLUMNS = [
    "project", "cancer_type", "status",
    "files_at_gdc", "files_downloaded", "files_incomplete", "md5_verified",
    "samples_total", "primary_tumor", "solid_tissue_normal", "metastatic",
    "patients_unique", "excluded_total", "excluded_detail",
    "genes", "cldn1_detected", "cldn1_median", "cldn1_min", "cldn1_max",
    "meta_aligned", "clinical", "size_gb",
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

    print("  Reading GTEx...")
    gtex = collect_gtex()
    normals = normals_table(df, gtex)
    print("  Reading GEO...")
    geo, geo_detail, geo_decisions = collect_geo()

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    root = os.path.dirname(DATASET_DIR)

    csv_path = os.path.join(SUMMARY_DIR, "tcga_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n  {os.path.relpath(csv_path, root)}")

    if not gtex.empty:
        gtex_csv = os.path.join(SUMMARY_DIR, "gtex_summary.csv")
        gtex.to_csv(gtex_csv, index=False)
        print(f"  {os.path.relpath(gtex_csv, root)}")

    if not geo.empty:
        geo_csv = os.path.join(SUMMARY_DIR, "geo_summary.csv")
        geo.to_csv(geo_csv, index=False)
        print(f"  {os.path.relpath(geo_csv, root)}")

    xlsx_path = os.path.join(SUMMARY_DIR, "wp1_summary.xlsx")
    if write_excel(xlsx_path, {"TCGA": df, "GTEx": gtex,
                               "Normals per cancer": normals,
                               "GEO candidates": geo,
                               "GEO inspected": geo_detail,
                               "GEO decisions": geo_decisions}):
        print(f"  {os.path.relpath(xlsx_path, root)}")
    else:
        print("  Excel skipped (pip install openpyxl)")

    def show(frame, cols):
        view = frame[[c for c in cols if c in frame.columns]]
        return view.apply(lambda c: ["-" if pd.isna(v) else str(v) for v in c])

    print("\n  === TCGA ===")
    print(show(df, ["project", "status", "samples_total", "primary_tumor",
                    "solid_tissue_normal", "metastatic", "excluded_total",
                    "cldn1_median"]).to_string(index=False))

    if not gtex.empty:
        print("\n  === GTEx (normal tissue) ===")
        print(show(gtex, ["tissue", "wp1_label", "samples", "donors",
                          "cldn1_median"]).to_string(index=False))

    if not geo.empty:
        print("\n  === GEO candidates ===")
        print(show(geo, ["cancer", "candidates_found", "past_gates_1_2",
                         "with_raw_counts", "tissue_confirmed", "ambiguous",
                         "downloaded"]).to_string(index=False))

    if not normals.empty:
        print("\n  === Normals per cancer (why GTEx is needed) ===")
        print(show(normals, ["project", "tcga_tumour", "tcga_normal",
                             "gtex_normal", "normal_total"]).to_string(index=False))


if __name__ == "__main__":
    start_logging()
    main()
