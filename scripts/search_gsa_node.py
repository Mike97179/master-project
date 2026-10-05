"""
search_gsa_node.py — CLAUD-IA WP1: Search Chinese genomic databases.

Usage:
    python scripts/search_gsa_node.py PAAD
    python scripts/search_gsa_node.py all

Searches:
    1. GSA (Genome Sequence Archive) — https://ngdc.cncb.ac.cn/gsa/
    2. NODE (National Omics Data Encyclopedia) — https://www.biosino.org/node/

Per Abrar's directive, these are secondary sources in case TCGA + GEO
don't provide enough samples.

Note: GSA has a public API; NODE search is more manual. This script
queries what it can and produces a CSV of candidate datasets to review.
"""

import os
import sys
import time
import requests
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES
from run_logger import start_logging

GSA_API = "https://ngdc.cncb.ac.cn/gsa/browse"


def search_gsa(keywords, max_results=50):
    """Search GSA for RNA-seq datasets. GSA uses a web API."""
    results = []

    for kw in keywords:
        try:
            r = requests.get(
                "https://ngdc.cncb.ac.cn/gsa/search",
                params={"searchTerm": f"{kw} RNA-seq homo sapiens", "pageSize": max_results},
                timeout=30,
                headers={"Accept": "application/json"},
            )
            if r.status_code == 200:
                try:
                    data = r.json()
                    if isinstance(data, dict) and "data" in data:
                        for item in data["data"]:
                            results.append({
                                "accession": item.get("accession", ""),
                                "title": item.get("title", "")[:120],
                                "organism": item.get("organism", ""),
                                "samples": item.get("sampleCount", "?"),
                                "repository": "GSA",
                            })
                except Exception:
                    pass
        except requests.exceptions.RequestException as e:
            print(f"    GSA request failed for '{kw}': {e}")

        time.sleep(0.5)

    return results


def search_node(keywords):
    """Search NODE. NODE's API is less documented, try basic search."""
    results = []

    for kw in keywords:
        try:
            r = requests.get(
                "https://www.biosino.org/node/api/search",
                params={"q": f"{kw} RNA-seq human", "type": "project", "size": 20},
                timeout=30,
                headers={"Accept": "application/json"},
            )
            if r.status_code == 200:
                try:
                    data = r.json()
                    items = data if isinstance(data, list) else data.get("hits", data.get("data", []))
                    for item in items:
                        if isinstance(item, dict):
                            results.append({
                                "accession": item.get("accession", item.get("id", "")),
                                "title": str(item.get("title", item.get("name", "")))[:120],
                                "organism": item.get("organism", ""),
                                "samples": item.get("sampleCount", "?"),
                                "repository": "NODE",
                            })
                except Exception:
                    pass
        except requests.exceptions.RequestException as e:
            print(f"    NODE request failed for '{kw}': {e}")

        time.sleep(0.5)

    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/search_gsa_node.py <CANCER|all>")
        sys.exit(1)

    target = sys.argv[1].upper()

    if target == "ALL":
        cancers = CANCER_TYPES
    else:
        cancers = {k: v for k, v in CANCER_TYPES.items() if v["short"] == target}
        if not cancers:
            print(f"Unknown: {target}. Valid: {[v['short'] for v in CANCER_TYPES.values()]}")
            sys.exit(1)

    all_results = []

    for project_id, info in cancers.items():
        print(f"\nSearching for: {info['label']}...")

        print("  Querying GSA...")
        gsa = search_gsa(info["keywords"])
        for r in gsa:
            r["cancer"] = info["short"]
        all_results.extend(gsa)
        print(f"  GSA: {len(gsa)} results")

        print("  Querying NODE...")
        node = search_node(info["keywords"])
        for r in node:
            r["cancer"] = info["short"]
        all_results.extend(node)
        print(f"  NODE: {len(node)} results")

        time.sleep(1)

    if all_results:
        df = pd.DataFrame(all_results)
        # Filter to human only
        if "organism" in df.columns:
            human_mask = df["organism"].str.contains("sapiens|human|homo", case=False, na=True)
            non_human = (~human_mask).sum()
            if non_human > 0:
                print(f"\n  Filtered out {non_human} non-human datasets")
            df = df[human_mask]

        df = df.drop_duplicates(subset="accession")

        print(f"\n{'='*80}")
        print(f"  Found {len(df)} datasets from Chinese databases")
        print(f"{'='*80}\n")

        if not df.empty:
            print(df[["repository", "accession", "cancer", "samples", "title"]].to_string(index=False))

        out_path = os.path.join(DATASET_DIR, "gsa_node_results.csv")
        os.makedirs(DATASET_DIR, exist_ok=True)
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path}")
    else:
        print("\nNo results found. These databases may need manual search:")
        print("  GSA:  https://ngdc.cncb.ac.cn/gsa/")
        print("  NODE: https://www.biosino.org/node/")


if __name__ == "__main__":
    start_logging()
    main()
