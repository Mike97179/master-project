"""
download_tcga.py — CLAUD-IA WP1: Download bulk RNA-seq from TCGA/GDC.

Usage:
    python scripts/download_tcga.py TCGA-PAAD
    python scripts/download_tcga.py all
    python scripts/download_tcga.py TCGA-PAAD --verify

Downloads STAR-Counts RNA-seq + clinical metadata for one or all 6 GI
cancer types. Requires gdc-client in PATH or in tools/.

Sample filter (Abrar directive):
    - Human biopsy only (Primary Tumor 01, Solid Tissue Normal 11, Metastatic 06)
    - No cell lines, no organoids

Output:
    dataset/TCGA/<PROJECT>/
    ├── rna/
    │   ├── manifest.txt
    │   ├── metadata.json
    │   └── raw/  (one folder per file UUID)
    ├── clinical/
    │   ├── clinical.tsv
    │   ├── follow_up.tsv
    │   ├── pathology_detail.tsv
    │   └── exposure.tsv
    └── processed/
        ├── counts_matrix.csv
        └── sample_metadata.csv
"""

import os
import sys
import json
import subprocess
import requests
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import (
    DATASET_DIR, CANCER_TYPES, GDC_FILES_URL, GDC_CASES_URL,
    TCGA_RNASEQ_FILTER, TCGA_VALID_SAMPLE_CODES, TCGA_SAMPLE_LABELS,
)
from run_logger import start_logging

# ── Arguments ─────────────────────────────────────────────────────────────

def parse_args():
    if len(sys.argv) < 2:
        print("Usage: python scripts/download_tcga.py <TCGA-XXXX|all> [--verify]")
        sys.exit(1)

    target = sys.argv[1]
    verify = "--verify" in sys.argv

    if target.lower() == "all":
        projects = list(CANCER_TYPES.keys())
    else:
        target = target.upper()
        if target not in CANCER_TYPES:
            print(f"ERROR: {target} is not a WP1 cancer type.")
            print(f"Valid: {', '.join(CANCER_TYPES.keys())}")
            sys.exit(1)
        projects = [target]

    return projects, verify

# ── GDC API helpers ───────────────────────────────────────────────────────

def build_gdc_filter(project_id):
    f = TCGA_RNASEQ_FILTER
    return {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_category", "value": f["data_category"]}},
            {"op": "=", "content": {"field": "data_type", "value": f["data_type"]}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": f["workflow_type"]}},
            {"op": "=", "content": {"field": "access", "value": f["access"]}},
        ],
    }


def query_gdc_files(project_id):
    filters = build_gdc_filter(project_id)
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,data_type,cases.submitter_id,"
                  "associated_entities.entity_submitter_id",
        "format": "json",
        "size": "2000",
    }
    r = requests.get(GDC_FILES_URL, params=params)
    r.raise_for_status()
    return r.json()["data"]


def download_manifest(project_id, manifest_path):
    filters = build_gdc_filter(project_id)
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,md5sum,file_size,state",
        "format": "json",
        "size": "2000",
    }
    r = requests.get(GDC_FILES_URL, params=params)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]

    lines = ["id\tfilename\tmd5\tsize\tstate"]
    for h in hits:
        lines.append(f"{h['file_id']}\t{h['file_name']}\t"
                      f"{h.get('md5sum','')}\t{h.get('file_size','')}\t"
                      f"{h.get('state','released')}")

    with open(manifest_path, "w") as f:
        f.write("\n".join(lines))
    return len(hits)


def download_clinical(project_id, clinical_dir):
    os.makedirs(clinical_dir, exist_ok=True)
    filters = {"op": "=", "content": {"field": "project.project_id", "value": project_id}}

    clinical_types = {
        "clinical":         "diagnoses,demographic",
        "follow_up":        "follow_ups",
        "pathology_detail": "diagnoses.pathology_details",
        "exposure":         "exposures",
    }

    for ctype, expand in clinical_types.items():
        out_path = os.path.join(clinical_dir, f"{ctype}.tsv")
        if os.path.exists(out_path):
            print(f"    {ctype}.tsv exists — skip")
            continue

        r = requests.get(GDC_CASES_URL, params={
            "filters": json.dumps(filters),
            "format": "tsv",
            "expand": expand,
            "size": "2000",
        })
        r.raise_for_status()
        with open(out_path, "w") as f:
            f.write(r.text)
        print(f"    {ctype}.tsv saved")


def run_gdc_client(manifest_path, output_dir):
    gdc_bin = "gdc-client"
    local_bin = os.path.join(os.path.dirname(__file__), "..", "tools", "gdc-client")
    if os.path.isfile(local_bin):
        gdc_bin = os.path.abspath(local_bin)

    cmd = [gdc_bin, "download", "-m", manifest_path, "-d", output_dir, "--n-processes", "4"]
    print(f"    $ {' '.join(cmd)}")

    # Stream through sys.stdout (not the inherited fd) so the run log
    # captures gdc-client's output as well as our own prints.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        data = os.read(proc.stdout.fileno(), 4096)
        if not data:
            break
        sys.stdout.write(data.decode("utf-8", "replace"))
    proc.stdout.close()
    return proc.wait()

# ── Assembly: raw TSV files → counts matrix ──────────────────────────────

def build_barcode_map(metadata):
    """Map file_name → TCGA barcode, filtering to valid sample types."""
    fmap = {}
    skipped = 0
    for record in metadata:
        barcode = record["associated_entities"][0]["entity_submitter_id"]
        sample_code = barcode.split("-")[3][:2]
        if sample_code not in TCGA_VALID_SAMPLE_CODES:
            skipped += 1
            continue
        fmap[record["file_name"]] = {
            "barcode": barcode,
            "patient": barcode[:12],
            "sample_type_code": sample_code,
            "sample_type_label": TCGA_SAMPLE_LABELS.get(sample_code, "Other"),
        }
    if skipped:
        print(f"    Filtered out {skipped} non-biopsy samples (cell lines, etc.)")
    return fmap


def assemble_counts(raw_dir, fname_map, out_dir):
    """Read individual STAR-Counts TSVs → single counts matrix."""
    tsv_files = []
    for root, _, files in os.walk(raw_dir):
        for f in files:
            if f.endswith(".tsv"):
                tsv_files.append(os.path.join(root, f))

    if not tsv_files:
        print("    No TSV files found — skipping assembly")
        return None, None

    counts = {}
    for filepath in sorted(tsv_files):
        fname = os.path.basename(filepath)
        if fname not in fname_map:
            continue
        barcode = fname_map[fname]["barcode"]
        df = pd.read_csv(filepath, sep="\t", comment="#",
                         usecols=["gene_id", "gene_name", "gene_type", "unstranded"])
        df = df[~df["gene_id"].str.startswith("N_")]
        df = df.set_index("gene_id")
        counts[barcode] = df["unstranded"]

    if not counts:
        print("    No samples matched metadata — skipping")
        return None, None

    matrix = pd.DataFrame(counts)
    ref = pd.read_csv(tsv_files[0], sep="\t", comment="#",
                      usecols=["gene_id", "gene_name", "gene_type"])
    ref = ref[~ref["gene_id"].str.startswith("N_")].set_index("gene_id")
    matrix.insert(0, "gene_name", ref["gene_name"])
    matrix.insert(1, "gene_type", ref["gene_type"])

    os.makedirs(out_dir, exist_ok=True)
    matrix.to_csv(os.path.join(out_dir, "counts_matrix.csv"))
    print(f"    Counts matrix: {matrix.shape[0]} genes x {matrix.shape[1]-2} samples")

    # Sample metadata
    sample_meta = pd.DataFrame([
        v for v in fname_map.values() if v["barcode"] in counts
    ])

    clinical_path = os.path.join(os.path.dirname(out_dir), "clinical", "clinical.tsv")
    if os.path.exists(clinical_path):
        clin = pd.read_csv(clinical_path, sep="\t", low_memory=False)
        clin = clin.drop_duplicates(subset="submitter_id")
        merge_cols = [c for c in [
            "submitter_id",
            "diagnoses.0.primary_diagnosis",
            "diagnoses.0.ajcc_pathologic_stage",
            "demographic.vital_status",
            "demographic.days_to_death",
            "demographic.age_at_index",
            "demographic.gender",
        ] if c in clin.columns]
        sample_meta = sample_meta.merge(
            clin[merge_cols], left_on="patient", right_on="submitter_id", how="left"
        ).drop(columns="submitter_id", errors="ignore")

    sample_meta.to_csv(os.path.join(out_dir, "sample_metadata.csv"), index=False)
    print(f"    Sample metadata: {len(sample_meta)} rows")

    return matrix, sample_meta

# ── Sense check ──────────────────────────────────────────────────────────

def sense_check(matrix, project_id):
    if matrix is None:
        return
    cldn1 = matrix[matrix["gene_name"] == "CLDN1"]
    if cldn1.empty:
        print("    WARNING: CLDN1 not found")
        return
    vals = cldn1.drop(columns=["gene_name", "gene_type"]).iloc[0]
    print(f"    CLDN1: detected in {(vals > 0).sum()}/{len(vals)} samples, "
          f"median={vals.median():.0f}, range=[{vals.min():.0f}, {vals.max():.0f}]")

# ── Verify mode ──────────────────────────────────────────────────────────

def verify_project(project_id):
    base = os.path.join(DATASET_DIR, "TCGA", project_id)
    data = query_gdc_files(project_id)
    expected = data["pagination"]["total"]
    raw_dir = os.path.join(base, "rna", "raw")
    downloaded = len([d for d in os.listdir(raw_dir) if os.path.isdir(os.path.join(raw_dir, d))]) if os.path.isdir(raw_dir) else 0
    processed = os.path.exists(os.path.join(base, "processed", "counts_matrix.csv"))
    status = "OK" if downloaded >= expected and processed else "INCOMPLETE"
    print(f"  {project_id}: API={expected}, downloaded={downloaded}, "
          f"processed={'yes' if processed else 'no'} → {status}")

# ── Main ─────────────────────────────────────────────────────────────────

def process_project(project_id):
    label = CANCER_TYPES[project_id]["label"]
    print(f"\n{'='*60}")
    print(f"  {project_id} — {label}")
    print(f"{'='*60}")

    base = os.path.join(DATASET_DIR, "TCGA", project_id)
    rna_dir = os.path.join(base, "rna")
    raw_dir = os.path.join(rna_dir, "raw")
    clinical_dir = os.path.join(base, "clinical")
    processed_dir = os.path.join(base, "processed")
    manifest_path = os.path.join(rna_dir, "manifest.txt")
    metadata_path = os.path.join(rna_dir, "metadata.json")

    os.makedirs(raw_dir, exist_ok=True)

    # 1. Query GDC
    print("\n  [1/5] Querying GDC API...")
    data = query_gdc_files(project_id)
    file_count = data["pagination"]["total"]
    print(f"    Files found: {file_count}")

    with open(metadata_path, "w") as f:
        json.dump(data["hits"], f, indent=2)

    # 2. Clinical
    print("\n  [2/5] Clinical metadata...")
    download_clinical(project_id, clinical_dir)

    # 3. Manifest + download
    print("\n  [3/5] Downloading RNA-seq files...")
    # Count only folders with a complete .tsv (not .partial)
    complete = 0
    for d in os.listdir(raw_dir):
        dpath = os.path.join(raw_dir, d)
        if not os.path.isdir(dpath):
            continue
        has_complete = any(f.endswith(".tsv") and not f.endswith(".partial")
                          for f in os.listdir(dpath))
        if has_complete:
            complete += 1
    if complete >= file_count:
        print(f"    Already downloaded ({complete} complete files) — skip")
    else:
        print(f"    Complete: {complete}/{file_count} — downloading missing...")
        n = download_manifest(project_id, manifest_path)
        print(f"    Manifest: {n} files")
        run_gdc_client(manifest_path, raw_dir)

    # 4. Assemble
    print("\n  [4/5] Assembling counts matrix...")
    fname_map = build_barcode_map(data["hits"])
    print(f"    Valid biopsy samples: {len(fname_map)}")
    matrix, sample_meta = assemble_counts(raw_dir, fname_map, processed_dir)

    # 5. Sense check
    print("\n  [5/5] Sense check...")
    sense_check(matrix, project_id)

    print(f"\n  Done: {base}")


if __name__ == "__main__":
    projects, verify = parse_args()
    start_logging()

    if verify:
        print("Verifying TCGA downloads...")
        for p in projects:
            verify_project(p)
    else:
        for p in projects:
            process_project(p)

    print("\n" + "="*60)
    print("  All TCGA downloads complete.")
    print("="*60)
