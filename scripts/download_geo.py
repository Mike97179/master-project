"""
download_geo.py — CLAUD-IA WP1: Download bulk RNA-seq datasets from GEO.

Usage:
    python scripts/download_geo.py GSE123456
    python scripts/download_geo.py search PAAD        # search for datasets
    python scripts/download_geo.py search all         # search all 6 cancers

Two modes:
    1. download <GSE_ID>  — download a specific dataset
    2. search <cancer|all> — search GEO for bulk RNA-seq datasets

Abrar filters applied:
    - Platform: RNA-seq only (GPL11154, GPL16791, GPL18573, GPL20301, etc.)
    - Organism: Homo sapiens
    - Sample type: tissue biopsy (flags cell lines and organoids in output)

Output:
    dataset/GEO/<GSE_ID>_<CANCER>/
    ├── expression.tsv
    ├── sample_metadata.tsv
    ├── summary.md
    └── raw/
"""

import os
import sys
import re
import shutil
import subprocess
import time
import requests
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES
from run_logger import start_logging

GEO_ROOT = os.path.join(DATASET_DIR, "GEO")

# ── Cancer type inference ────────────────────────────────────────────────
KEYWORD_TO_CANCER = {}
for project_id, info in CANCER_TYPES.items():
    short = info["short"]
    for kw in info["keywords"]:
        KEYWORD_TO_CANCER[kw.lower()] = short

KEYWORD_TO_CANCER.update({
    "cell line encyclopedia": "CELL_LINE",
    "ccle": "CELL_LINE",
    "organoid": "ORGANOID",
})

# ── Excluded sample patterns (Abrar: no cell lines, no organoids) ────────
EXCLUDED_PATTERNS = re.compile(
    r"cell\s*line|organoid|xenograft|pdx|mouse|murine|rattus|"
    r"mus musculus|drosophila|zebrafish|caenorhabditis",
    re.IGNORECASE,
)


def infer_cancer(title, summary):
    text = (title + " " + summary).lower()
    matches = set()
    for kw, cancer in KEYWORD_TO_CANCER.items():
        if kw in text:
            matches.add(cancer)
    matches.discard("CELL_LINE")
    matches.discard("ORGANOID")
    if len(matches) == 1:
        return matches.pop()
    if not matches:
        return "UNKNOWN"
    for kw, cancer in KEYWORD_TO_CANCER.items():
        if kw in title.lower() and cancer not in ("CELL_LINE", "ORGANOID"):
            return cancer
    return "AMBIGUOUS"


def flag_sample_issues(title, summary):
    """Check if dataset likely contains non-biopsy samples."""
    text = (title + " " + summary).lower()
    flags = []
    if "cell line" in text:
        flags.append("CELL_LINE")
    if "organoid" in text:
        flags.append("ORGANOID")
    if "xenograft" in text or "pdx" in text:
        flags.append("XENOGRAFT")
    if any(w in text for w in ["mouse", "murine", "mus musculus"]):
        flags.append("NON_HUMAN")
    return flags


# ── NCBI E-utilities ─────────────────────────────────────────────────────

def fetch_geo_meta(gse_id):
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                     params={"db": "gds", "term": f"{gse_id}[ACCN]", "retmode": "json"})
    ids = r.json()["esearchresult"].get("idlist", [])
    if not ids:
        raise ValueError(f"No NCBI record for {gse_id}")
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                     params={"db": "gds", "id": ids[0], "retmode": "json"})
    return r.json()["result"].get(ids[0], {})


def search_geo_rnaseq(keywords, max_results=50):
    """Search GEO for bulk RNA-seq datasets matching keywords."""
    query = (f'({" OR ".join(keywords)}) AND '
             '"expression profiling by high throughput sequencing"[DataSet Type] AND '
             '"homo sapiens"[Organism]')

    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                     params={"db": "gds", "term": query, "retmax": max_results,
                             "retmode": "json"})
    ids = r.json()["esearchresult"].get("idlist", [])
    if not ids:
        return []

    results = []
    for gds_id in ids:
        time.sleep(0.35)
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                         params={"db": "gds", "id": gds_id, "retmode": "json"})
        meta = r.json()["result"].get(gds_id, {})
        if not meta:
            continue
        accession = meta.get("accession", "")
        if not accession.startswith("GSE"):
            continue

        title = meta.get("title", "")
        summary = meta.get("summary", "")
        flags = flag_sample_issues(title, summary)

        results.append({
            "accession": accession,
            "title": title[:100],
            "n_samples": meta.get("n_samples", "?"),
            "cancer": infer_cancer(title, summary),
            "flags": ", ".join(flags) if flags else "OK",
            "gpl": meta.get("gpl", "?"),
            "date": meta.get("pdat", "?"),
        })

    return results


# ── Download a single GSE ────────────────────────────────────────────────

def download_gse(gse_id):
    print(f"\n{'='*60}")
    print(f"  Downloading {gse_id}")
    print(f"{'='*60}")

    # 1. Metadata
    print("\n  [1/4] Querying NCBI...")
    meta = fetch_geo_meta(gse_id)
    title = meta.get("title", "")
    summary = meta.get("summary", "")
    cancer = infer_cancer(title, summary)
    flags = flag_sample_issues(title, summary)

    print(f"    Title:  {title[:80]}")
    print(f"    Cancer: {cancer}")
    if flags:
        print(f"    WARNING: sample issues detected: {', '.join(flags)}")

    # 2. Create folder
    folder = f"{gse_id}_{cancer}"
    out_dir = os.path.join(GEO_ROOT, folder)
    raw_dir = os.path.join(out_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    # 3. Download via GEOparse
    print("\n  [2/4] Downloading series matrix...")
    try:
        import GEOparse
        tmp_dir = os.path.join(out_dir, "_tmp")
        os.makedirs(tmp_dir, exist_ok=True)
        gse = GEOparse.get_GEO(geo=gse_id, destdir=tmp_dir, silent=True)

        # Expression matrix
        print("  [3/4] Extracting expression...")
        expr_dfs = []
        for gsm_id, gsm in gse.gsms.items():
            if gsm.table is None or gsm.table.empty:
                continue
            df = gsm.table.set_index(gsm.table.columns[0])
            expr_dfs.append(df[df.columns[0]].rename(gsm_id))

        expr_path = os.path.join(out_dir, "expression.tsv")
        if expr_dfs:
            expr = pd.concat(expr_dfs, axis=1)
            expr.to_csv(expr_path, sep="\t")
            print(f"    Expression: {expr.shape[0]} features x {expr.shape[1]} samples")
        else:
            expr = None
            print("    WARNING: no expression data in GSMs")

        # Sample metadata
        rows = []
        for gsm_id, gsm in gse.gsms.items():
            row = {"gsm_id": gsm_id}
            for key, vals in gsm.metadata.items():
                row[key] = "; ".join(str(v) for v in vals) if isinstance(vals, list) else str(vals)
            rows.append(row)
        sample_df = pd.DataFrame(rows)
        sample_df.to_csv(os.path.join(out_dir, "sample_metadata.tsv"), sep="\t", index=False)
        print(f"    Sample metadata: {len(sample_df)} samples")

        shutil.rmtree(tmp_dir, ignore_errors=True)
    except Exception as e:
        print(f"    ERROR: {e}")
        expr = None

    # 4. Download raw supplementary files
    print("\n  [4/4] Downloading raw supplementary files...")
    stub = gse_id[:-3] + "nnn"
    ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse_id}/suppl/"
    cmd = ["wget", "--no-verbose", "--no-host-directories", "--cut-dirs=5",
           "--recursive", "--no-parent", "--reject", "index.html*,robots.txt",
           "--directory-prefix", raw_dir, ftp_url]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"    wget failed — download manually from:")
        print(f"    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}")

    # Summary
    summary_path = os.path.join(out_dir, "summary.md")
    with open(summary_path, "w") as f:
        f.write(f"# {gse_id} — {title}\n\n")
        f.write(f"- **Cancer:** {cancer}\n")
        f.write(f"- **Samples:** {meta.get('n_samples', '?')}\n")
        f.write(f"- **Platform:** GPL{meta.get('gpl', '?')}\n")
        f.write(f"- **Sample issues:** {', '.join(flags) if flags else 'None (human biopsy)'}\n")
        f.write(f"- **Date:** {meta.get('pdat', '?')}\n\n")
        f.write(f"## Description\n{meta.get('summary', '')}\n\n")
        f.write(f"## Source\nhttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}\n")

    print(f"\n  Done: {out_dir}")


# ── Main ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python scripts/download_geo.py GSE123456")
        print("  python scripts/download_geo.py search PAAD")
        print("  python scripts/download_geo.py search all")
        sys.exit(1)

    start_logging()
    mode = sys.argv[1].lower()

    if mode == "search":
        target = sys.argv[2].upper() if len(sys.argv) > 2 else "ALL"

        if target == "ALL":
            cancers = CANCER_TYPES
        else:
            cancers = {k: v for k, v in CANCER_TYPES.items() if v["short"] == target}
            if not cancers:
                print(f"Unknown cancer: {target}. Valid: {[v['short'] for v in CANCER_TYPES.values()]}")
                sys.exit(1)

        all_results = []
        for project_id, info in cancers.items():
            print(f"\nSearching GEO for: {info['label']}...")
            results = search_geo_rnaseq(info["keywords"])
            for r in results:
                r["search_cancer"] = info["short"]
            all_results.extend(results)
            time.sleep(1)

        if all_results:
            df = pd.DataFrame(all_results)
            df = df.drop_duplicates(subset="accession")
            print(f"\n{'='*80}")
            print(f"  Found {len(df)} unique datasets")
            print(f"{'='*80}\n")
            print(df[["accession", "cancer", "n_samples", "flags", "title"]].to_string(index=False))

            out_path = os.path.join(GEO_ROOT, "search_results.csv")
            os.makedirs(GEO_ROOT, exist_ok=True)
            df.to_csv(out_path, index=False)
            print(f"\nSaved: {out_path}")
        else:
            print("No results found.")

    elif mode.startswith("gse"):
        download_gse(mode.upper())

    else:
        download_gse(mode.upper())
