"""
search_pubmed.py — CLAUD-IA WP1: Find bulk RNA-seq datasets via PubMed.

Usage:
    python scripts/search_pubmed.py PAAD
    python scripts/search_pubmed.py all
    python scripts/search_pubmed.py all --full

Searches PubMed for papers that mention RNA-seq + a GI cancer type, then
checks whether the abstract or MeSH terms mention data repositories
(GEO, ArrayExpress, SRA, GSA, NODE, Zenodo, etc.).

Per Abrar's directive: "in the 'data availability' section, they put
where they have the data."

Output:
    dataset/pubmed_datasets.csv — table of papers with data pointers
"""

import os
import sys
import time
import re
import requests
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES
from run_logger import start_logging

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

REPO_PATTERNS = [
    (r"GSE\d{4,8}", "GEO"),
    (r"PRJNA\d{5,9}", "SRA/BioProject"),
    (r"SRP\d{5,9}", "SRA"),
    (r"E-MTAB-\d{4,8}", "ArrayExpress"),
    (r"CRA\d{5,9}", "GSA (China)"),
    (r"PRJCA\d{5,9}", "GSA (China)"),
    (r"OEP\d{5,9}", "NODE (China)"),
    (r"PRJEB\d{5,9}", "ENA"),
    (r"10\.\d{4,}/zenodo\.\d+", "Zenodo"),
    (r"synapse\.org", "Synapse"),
]


def search_pubmed(keywords, max_results=100):
    query = (f'({" OR ".join(keywords)}) AND '
             '(RNA-seq OR "RNA sequencing" OR "bulk RNA") AND '
             '(dataset OR "data availability" OR "deposited" OR "accession")')

    r = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params={
        "db": "pubmed", "term": query, "retmax": max_results,
        "retmode": "json", "sort": "relevance",
    })
    ids = r.json()["esearchresult"].get("idlist", [])
    return ids


def fetch_abstracts(pmids, batch_size=50):
    results = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        r = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params={
            "db": "pubmed", "id": ",".join(batch), "rettype": "abstract",
            "retmode": "xml",
        })
        results.append(r.text)
        time.sleep(0.4)
    return "\n".join(results)


def extract_repos_from_text(text):
    found = []
    for pattern, repo_name in REPO_PATTERNS:
        matches = re.findall(pattern, text, re.IGNORECASE)
        for m in matches:
            found.append({"accession": m, "repository": repo_name})
    return found


def parse_articles(xml_text):
    """Quick regex-based parse to avoid lxml dependency."""
    articles = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml_text, re.DOTALL)
    rows = []
    for art in articles:
        pmid = re.search(r"<PMID.*?>(\d+)</PMID>", art)
        title = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", art, re.DOTALL)
        abstract = re.search(r"<AbstractText.*?>(.*?)</AbstractText>", art, re.DOTALL)
        year = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", art, re.DOTALL)

        pmid = pmid.group(1) if pmid else ""
        title_text = title.group(1) if title else ""
        abstract_text = abstract.group(1) if abstract else ""
        year_text = year.group(1) if year else ""

        full_text = title_text + " " + abstract_text + " " + art
        repos = extract_repos_from_text(full_text)

        if repos:
            for repo in repos:
                rows.append({
                    "pmid": pmid,
                    "year": year_text,
                    "title": re.sub(r"<.*?>", "", title_text)[:120],
                    "accession": repo["accession"],
                    "repository": repo["repository"],
                })

    return rows


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/search_pubmed.py <CANCER|all> [--full]")
        sys.exit(1)

    target = sys.argv[1].upper()
    full_mode = "--full" in sys.argv
    max_results = 200 if full_mode else 50

    if target == "ALL":
        cancers = CANCER_TYPES
    else:
        cancers = {k: v for k, v in CANCER_TYPES.items() if v["short"] == target}
        if not cancers:
            print(f"Unknown: {target}. Valid: {[v['short'] for v in CANCER_TYPES.values()]}")
            sys.exit(1)

    all_rows = []

    for project_id, info in cancers.items():
        print(f"\nSearching PubMed for: {info['label']}...")
        pmids = search_pubmed(info["keywords"], max_results=max_results)
        print(f"  Found {len(pmids)} papers")

        if pmids:
            xml = fetch_abstracts(pmids)
            rows = parse_articles(xml)
            for r in rows:
                r["cancer"] = info["short"]
            all_rows.extend(rows)
            print(f"  Papers with data accessions: {len(set(r['pmid'] for r in rows))}")

        time.sleep(1)

    if all_rows:
        df = pd.DataFrame(all_rows).drop_duplicates(subset=["accession", "cancer"])
        df = df.sort_values(["cancer", "repository", "accession"])

        print(f"\n{'='*80}")
        print(f"  Found {len(df)} dataset accessions across {df['cancer'].nunique()} cancer types")
        print(f"{'='*80}\n")

        summary = df.groupby(["cancer", "repository"]).size().reset_index(name="count")
        print(summary.to_string(index=False))

        print(f"\nTop accessions:")
        for _, row in df.head(20).iterrows():
            print(f"  {row['repository']:15s} {row['accession']:20s} {row['cancer']:6s} {row['title'][:60]}")

        out_path = os.path.join(DATASET_DIR, "pubmed_datasets.csv")
        os.makedirs(DATASET_DIR, exist_ok=True)
        df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path}")
    else:
        print("\nNo datasets found in abstracts.")


if __name__ == "__main__":
    start_logging()
    main()
