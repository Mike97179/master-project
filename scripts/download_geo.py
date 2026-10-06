"""
download_geo.py — CLAUD-IA WP1: find and download bulk RNA-seq from GEO.

Usage:
    python scripts/download_geo.py search PAAD      # search one cancer
    python scripts/download_geo.py search all       # search all 6
    python scripts/download_geo.py inspect GSE335154  # check before downloading
    python scripts/download_geo.py download GSE335154

Why three modes. TCGA and GTEx ship one uniform matrix each; GEO does not.
Every submitter uploads whatever they like, and many upload only normalised
values (CPM/FPKM/TPM), which DESeq2 cannot use — it needs raw counts. So the
workflow is search, then inspect the candidates, then download only the ones
that carry raw counts.

Abrar's filters:
    - Human tissue biopsy: cell lines, organoids, xenografts and non-human
      datasets are flagged, not silently dropped, because GEO's free-text
      metadata is not reliable enough to decide automatically
    - Bulk RNA-seq only (DataSet Type = expression profiling by high
      throughput sequencing)

Output:
    dataset/GEO/
    ├── inspected.csv                   verdict cache, shared across cancers
    └── <CANCER>/                       one folder per WP1 cancer type
        ├── search_results.csv          candidates this search found
        ├── decisions.csv               why each one was kept or skipped
        └── <GSE>/
            ├── summary.md
            ├── sample_metadata.tsv
            ├── suppl_files.csv
            └── raw/
"""

import os
import re
import shutil
import subprocess
import sys
import time

import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(__file__))
from config import DATASET_DIR, CANCER_TYPES
from run_logger import start_logging

GEO_ROOT = os.path.join(DATASET_DIR, "GEO")


def cancer_dir(short):
    """dataset/GEO/<CANCER>/ — each cancer's curation kept separate.

    The verdict cache stays at the GEO root instead, because a series'
    verdict does not depend on which cancer's keywords found it.
    """
    path = os.path.join(GEO_ROOT, short.upper())
    os.makedirs(path, exist_ok=True)
    return path
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# NCBI allows 3 requests/second without an API key
NCBI_PAUSE = 0.35

# ── Cancer type inference ────────────────────────────────────────────────
KEYWORD_TO_CANCER = {}
for project_id, info in CANCER_TYPES.items():
    for kw in info["keywords"]:
        KEYWORD_TO_CANCER[kw.lower()] = info["short"]

# ── Sample-type flags (Abrar: human biopsy only) ─────────────────────────
FLAG_PATTERNS = {
    "CELL_LINE": (r"cell\s*lines?|ccle|\bhek\b|\bhela\b|panc-?1|mia\s*paca|"
                  r"bxpc-?3|capan|aspc-?1|\bhct-?116\b|\bhepg2\b|\bsw480\b|"
                  r"\bshrna\b|\bsirna\b|knock-?down|knock-?out|\bko\b|"
                  r"[a-z0-9]{2,}ko\b|"
                  r"\bcrispr\b|transfect|overexpress|\bparental\b|\bclones?\b|"
                  # lines that recur in liver and biliary work
                  r"\blo2\b|\bthle-?2?\b|hucct|\brbe\b|qbc939|huh-?\d+|tfk-?1|"
                  r"\bccl?p-?1\b|snu-?\d{3,4}|egi-?1|mz-?cha|gbc-?sd|\bnoz\b|"
                  # shRNA/siRNA designs name their columns this way
                  r"\bsh\d\b|\bsi\d\b|\bscr\d?\b|scramble|\bshctrl\b|\bshnc\b|"
                  r"\bshluc\b|\bnegative control\b"),
    "ORGANOID": r"organoid|\bspheroid\b",
    "XENOGRAFT": r"xenograft|\bpdx\b",
    "NON_HUMAN": r"\bmouse\b|\bmice\b|murine|mus musculus|\brat\b|rattus|zebrafish|drosophila",
    "SINGLE_CELL": (r"single[-\s]?cell|single[-\s]?nucle|\bscrna\b|\bsnrna\b|"
                    r"10x genomics|fragments\.tsv|barcodes\.tsv|features\.tsv|matrix\.mtx|"
                    r"cell[-_\s]?annotation"),
    "NOT_RNASEQ": r"\batac\b|\bchip[-\s]?seq\b|\bcut&?run\b|bisulfite|\bwgbs\b|\bmethylat",
    # Liquid biopsy is not a tissue biopsy: GSE183635 ships 2351 samples of
    # platelet RNA, which passed every other filter.
    "NOT_TISSUE": (r"platelets?|liquid biops|\bplasma\b|\bserum\b|whole blood|"
                   r"\bpbmc\b|buffy coat|cell[-\s]?free|\bcfrna\b|\bctdna\b|"
                   r"exosom|extracellular vesicle|\bsaliva\b|\burine\b"),
}

# ── Supplementary file classification ────────────────────────────────────
# Checked in this order: a file called "Normcount" is normalised, not raw,
# so the normalised markers have to win over the word "count".
NORMALISED_MARKERS = r"fpkm|rpkm|\btpm\b|\bcpm\b|norm|scaled|vst|rlog|z[-_]?score|deseq|edger|\blog2\b"
RAW_MARKERS = r"raw|count|htseq|featurecount|star|\bexpected_count\b|quant"
# Counts that are not per gene are no use for a gene-level analysis:
# GSE138109's only "count" file is circRNAs_count, with gene data published
# as TPM alone.
NON_GENE_MARKERS = (r"\bcircrnas?\b|circular rna|\bisoforms?\b|\btranscripts?\b|"
                    r"\bmirnas?\b|microrna|\bexons?\b|\bjunctions?\b|\bpeaks?\b|"
                    r"\bpromoters?\b|\bsplic")


def classify_suppl(fname):
    """raw / normalised / bundle / other, from the file name alone.

    Underscores and dots become spaces first: regex treats "_" as a word
    character, so \bcpm\b would miss "CPM_Expression_Profile.xlsx".

    "<GSE>_RAW.tar" is a special case. It is GEO's generic container for
    every per-sample supplementary file, so the word RAW there means "the
    raw supplementary files", not "raw counts" — GSE348275_RAW.tar turned
    out to hold single-cell ATAC fragments. It gets its own class so it
    never counts as evidence of a usable count matrix.
    """
    low = fname.lower()
    if re.search(r"_raw\.(tar|zip)$", low):
        return "bundle"
    tokens = re.sub(r"[^a-z0-9]+", " ", low)
    if re.search(NON_GENE_MARKERS, tokens):
        return "non_gene"
    if re.search(NORMALISED_MARKERS, tokens):
        return "normalised"
    if re.search(RAW_MARKERS, tokens):
        return "raw"
    return "other"


def parse_size(text):
    """'1.8M' -> bytes. The HTTPS listing gives sizes in this short form."""
    m = re.match(r"([0-9.]+)\s*([KMGT]?)", (text or "").strip(), re.I)
    if not m:
        return 0
    mult = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}
    return int(float(m.group(1)) * mult[m.group(2).upper()])


def human_size(n):
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024


def infer_cancer(title, summary):
    text = (title + " " + summary).lower()
    hits = {c for kw, c in KEYWORD_TO_CANCER.items() if kw in text}
    if len(hits) == 1:
        return hits.pop()
    if not hits:
        return "UNKNOWN"
    # Several cancers mentioned: trust the title over the summary
    for kw, cancer in KEYWORD_TO_CANCER.items():
        if kw in title.lower():
            return cancer
    return "AMBIGUOUS"


def flag_sample_issues(title, summary):
    """Sample-type flags from free text or from a list of file names.

    The text is tokenised, because regex counts "_" as a word character:
    without this, \blo2\b never matches GSE276342_LO2_raw_count.txt.gz,
    and that series is the LO2 hepatocyte line.
    """
    raw = (title + " " + summary).lower()
    tokens = re.sub(r"[^a-z0-9&]+", " ", raw)
    return [name for name, pat in FLAG_PATTERNS.items()
            if re.search(pat, raw) or re.search(pat, tokens)]


# ── NCBI E-utilities ─────────────────────────────────────────────────────

def esummary(gds_id):
    r = requests.get(f"{EUTILS}/esummary.fcgi",
                     params={"db": "gds", "id": gds_id, "retmode": "json"},
                     timeout=60)
    r.raise_for_status()
    return r.json()["result"].get(gds_id, {})


def fetch_geo_meta(gse_id):
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gds", "term": f"{gse_id}[ACCN]", "retmode": "json"},
                     timeout=60)
    r.raise_for_status()
    ids = r.json()["esearchresult"].get("idlist", [])
    if not ids:
        raise ValueError(f"No NCBI record for {gse_id}")
    time.sleep(NCBI_PAUSE)
    return esummary(ids[0])


def search_geo_rnaseq(keywords):
    """Every candidate dataset for one cancer.

    There is no cap, on purpose. NCBI returns the most recent first, so a
    capped search is a biased sample rather than a smaller one: what gets
    published now is mostly single-cell and cell-line work, while bulk
    RNA-seq of biopsies is commoner among the older entries. A full sweep
    of cholangiocarcinoma found 22 usable datasets where the first 50
    entries held almost none.
    """
    query = (f'({" OR ".join(keywords)}) AND '
             '"expression profiling by high throughput sequencing"[DataSet Type] AND '
             '"homo sapiens"[Organism]')
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gds", "term": query, "retmax": 100000,
                             "retmode": "json"}, timeout=120)
    r.raise_for_status()
    res = r.json()["esearchresult"]
    total = int(res.get("count", 0))
    ids = res.get("idlist", [])
    print(f"    {total} datasets match; retrieving all {len(ids)}")
    mins = len(ids) * NCBI_PAUSE / 60
    if mins > 2:
        print(f"    ~{mins:.0f} min just to list them (NCBI allows 3 requests/second)")

    results = []
    for gds_id in ids:
        time.sleep(NCBI_PAUSE)
        meta = esummary(gds_id)
        accession = meta.get("accession", "")
        if not accession.startswith("GSE"):
            continue
        title = meta.get("title", "")
        summary = meta.get("summary", "")
        flags = flag_sample_issues(title, summary)
        results.append({
            "accession": accession,
            "inferred_cancer": infer_cancer(title, summary),
            "n_samples": meta.get("n_samples", ""),
            "flags": ", ".join(flags) if flags else "OK",
            "platform": f"GPL{meta.get('gpl', '?')}",
            "date": meta.get("pdat", ""),
            "title": title[:120],
            "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        })
    return results


# ── Supplementary file listing (no download) ─────────────────────────────

def list_suppl_files(gse_id):
    """Names and sizes of a series' supplementary files, over HTTPS.

    NCBI serves the FTP tree over HTTPS too, which is faster and less likely
    to be blocked than plain FTP.
    """
    stub = gse_id[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse_id}/suppl/"
    r = requests.get(url, timeout=90)
    if not r.ok:
        return [], url
    files = []
    # The listing is an HTML table: href, then date and size as plain text
    for m in re.finditer(r'href="([^"?/][^"]*)"[^<]*</a>\s*([0-9-]{10}\s[0-9:]{5})?\s*([0-9.]+[KMG]?)?',
                         r.text):
        name = m.group(1)
        if name.startswith("http") or name == "../":
            continue
        files.append({"file": name, "size": (m.group(3) or "").strip(),
                      "kind": classify_suppl(name)})
    return files, url


def peek_file(gse_id, fname, max_bytes=262144):
    """First few KB of a supplementary file, without downloading all of it.

    Uses an HTTP range request, decompressing on the fly for .gz. Needed
    because file names lie: a series shipping "Rawcount.csv.gz" full of
    decimals has not shipped raw counts.
    """
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/suppl/{fname}"
    try:
        r = requests.get(url, headers={"Range": f"bytes=0-{max_bytes}"},
                         timeout=90, stream=True)
        if r.status_code not in (200, 206):
            return ""
        data = r.content
        if fname.endswith(".gz"):
            import zlib
            d = zlib.decompressobj(16 + zlib.MAX_WBITS)
            data = d.decompress(data)
        return data.decode("utf-8", "replace")
    except Exception:
        return ""


def values_are_integers(text):
    """Are the numbers in this matrix whole? True / False / None if unknown.

    Tokenising beats a regex over the raw text. Fields can be separated by
    spaces as well as tabs or commas — GSE179443 ships a space-separated
    "Raw_gene_counts" whose values are decimals — and the first field of
    each row is a gene id whose version suffix (ENSG00000240361.2) looks
    exactly like a decimal number.

    A value of "5.0" still counts as whole: R writes integers that way.
    """
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    checked = 0
    for line in lines[1:60]:
        fields = [f.strip().strip('"').strip("'")
                  for f in re.split(r"[\t,;]|\s+", line.strip())]
        for v in fields[1:]:                       # field 0 is the gene id
            if not re.fullmatch(r"-?\d+(\.\d+)?([eE][-+]?\d+)?", v):
                continue
            checked += 1
            num = float(v)
            if num != int(num):
                return False
    return True if checked else None


def verify_contents(gse_id, files):
    """Look inside the files: are the 'raw' counts integers, and do the
    sample annotations betray cell lines the abstract never mentioned?

    Both filters upstream work on strings — file names and the free-text
    abstract — and both can be wrong in the same direction, towards a
    dataset looking more usable than it is.
    """
    out = {"raw_is_integer": "", "content_flags": ""}

    header_flags = []
    raw = next((f for f in files if f["kind"] == "raw"), None)
    if raw:
        text = peek_file(gse_id, raw["file"])
        if text:
            verdict = values_are_integers(text)
            out["raw_is_integer"] = {True: "yes", False: "no", None: ""}[verdict]
            # The column names of a count matrix name the samples, and they
            # are the frankest metadata of all: GSE312961's columns are
            # GBC-SD-CASE and NOZ-CASE, two gallbladder lines, and
            # GSE149536's are count.scr1..3 / count.sh1..3, a shRNA
            # knockdown. Neither abstract says "cell line".
            first_line = text.splitlines()[0] if text.splitlines() else ""
            header_flags = flag_sample_issues(first_line[:4000], "")

    # Small annotation files often name the cell lines outright, and a
    # bundle's filelist.txt names every file inside it, which is what
    # exposes single-cell or ATAC data hiding in a "_RAW.tar".
    flags = set()
    for f in files:
        if f["kind"] not in ("other", "bundle"):
            continue
        if f["kind"] == "bundle" or parse_size(f["size"]) > 2 * 1024 ** 2:
            continue
        text = peek_file(gse_id, f["file"], 65536).lower()
        for name, pat in FLAG_PATTERNS.items():
            if re.search(pat, text):
                flags.add(name)
    out["content_flags"] = ", ".join(sorted(set(flags) | set(header_flags)))
    return out


def inspect_gse(gse_id, quiet=False):
    """Report what a series contains, without downloading the data.

    quiet=True returns the verdict without the per-dataset report, so a
    batch run prints one line each instead of a screenful.
    """
    if not quiet:
        print(f"\n{'=' * 60}")
        print(f"  Inspecting {gse_id}")
        print(f"{'=' * 60}")

    meta = fetch_geo_meta(gse_id)
    title = meta.get("title", "")
    summary = meta.get("summary", "")
    cancer = infer_cancer(title, summary)
    flags = flag_sample_issues(title, summary)

    if not quiet:
        print(f"\n  Title:     {title[:90]}")
        print(f"  Samples:   {meta.get('n_samples', '?')}")
        print(f"  Platform:  GPL{meta.get('gpl', '?')}   Date: {meta.get('pdat', '?')}")
        print(f"  Cancer:    {cancer}")
        print(f"  Flags:     {', '.join(flags) if flags else 'OK (looks like human biopsy)'}")

    files, url = list_suppl_files(gse_id)

    # The file names are metadata too, and they are often more candid than
    # the abstract: GSE162739's abstract never says "cell line", but its
    # files are named Bap1KO-C91S-vs-PAR_diffexp, i.e. engineered cultures.
    name_flags = flag_sample_issues(" ".join(f["file"] for f in files), "")
    for fl in name_flags:
        if fl not in flags:
            flags.append(fl)

    if not quiet:
        print(f"\n  Supplementary files ({len(files)}):")
        if not files:
            print(f"    none listed — check {url}")
        for f in files:
            print(f"    [{f['kind']:10}] {f['size']:>7}  {f['file']}")

    kinds = {f["kind"] for f in files}
    if "raw" in kinds:
        verdict = "USABLE — raw counts available"
    elif "bundle" in kinds:
        verdict = "BUNDLE — per-sample archive, check filelist.txt by hand"
    elif "normalised" in kinds:
        verdict = "NORMALISED ONLY — no raw counts for DESeq2"
    elif "non_gene" in kinds:
        verdict = "NOT GENE-LEVEL — only circRNA/isoform/miRNA counts"
    elif files:
        verdict = "UNCLEAR — inspect the files by hand"
    else:
        verdict = "NO SUPPLEMENTARY FILES"
    if not quiet:
        print(f"\n  Verdict: {verdict}")
        if flags:
            print(f"  Note: flagged as {', '.join(flags)} — confirm sample types before use")

    row = {
        "accession": gse_id,
        "inferred_cancer": cancer,
        "n_samples": meta.get("n_samples", ""),
        "flags": ", ".join(flags) if flags else "OK",
        "platform": f"GPL{meta.get('gpl', '?')}",
        "n_suppl_files": len(files),
        "has_raw": "yes" if "raw" in kinds else "no",
        "has_normalised": "yes" if "normalised" in kinds else "no",
        "suppl_bytes": sum(parse_size(f["size"]) for f in files),
        # What a download would actually fetch: a _RAW.tar bundle is skipped
        # when the series also publishes a count matrix.
        "download_bytes": sum(
            parse_size(f["size"]) for f in files
            if not (f["kind"] == "bundle" and any(g["kind"] == "raw" for g in files))),
        "verdict": verdict,
        "title": title[:120],
        "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}",
    }
    record_inspection(row, quiet=quiet)
    return row, files


def record_inspection(row, quiet=False):
    """Accumulate verdicts in one CSV, newest value per accession."""
    os.makedirs(GEO_ROOT, exist_ok=True)
    path = os.path.join(GEO_ROOT, "inspected.csv")
    df = pd.DataFrame([row])
    if os.path.exists(path):
        old = pd.read_csv(path)
        df = pd.concat([old[old["accession"] != row["accession"]], df], ignore_index=True)
    df.to_csv(path, index=False)
    if not quiet:
        print(f"  Recorded in {os.path.relpath(path, os.path.dirname(DATASET_DIR))}")


def inspect_batch(short, ok_only=False, recheck=False):
    """Inspect every candidate from the last search, resuming where it left off.

    Inspecting by hand does not scale: a search returns 50 or more candidates
    and each one needs two NCBI calls. This walks the list, skips what is
    already in inspected.csv, and prints one line per dataset.
    """
    search_csv = os.path.join(cancer_dir(short), "search_results.csv")
    if not os.path.exists(search_csv):
        print(f"  No search results for {short} — run a search first:")
        print(f"    python scripts/download_geo.py search {short}")
        sys.exit(1)

    df = pd.read_csv(search_csv)
    if ok_only:
        df = df[df["flags"] == "OK"]
        print(f"  Restricted to the {len(df)} candidates flagged OK")

    insp_csv = os.path.join(GEO_ROOT, "inspected.csv")
    done = set()
    if os.path.exists(insp_csv) and not recheck:
        done = set(pd.read_csv(insp_csv)["accession"])

    todo = [a for a in df["accession"] if a not in done]
    print(f"  {len(df)} candidates, {len(df) - len(todo)} already inspected, "
          f"{len(todo)} to go")
    if not todo:
        print("  Nothing to do.")
        return

    verdicts = []
    for i, acc in enumerate(todo, start=1):
        time.sleep(NCBI_PAUSE)
        try:
            row, _ = inspect_gse(acc, quiet=True)
        except Exception as e:
            print(f"  [{i}/{len(todo)}] {acc}: ERROR {e}")
            continue
        verdicts.append(row)
        print(f"  [{i}/{len(todo)}] {acc:12} raw={row['has_raw']:3} "
              f"{row['flags'][:22]:22} {row['verdict']}")

    if not verdicts:
        return
    out = pd.DataFrame(verdicts)
    print(f"\n  === Verdicts ===")
    for v, n in out["verdict"].value_counts().items():
        print(f"    {n:4d}  {v}")
    usable = out[(out["has_raw"] == "yes") & (out["flags"] == "OK")]
    print(f"\n  Usable (raw counts AND no sample-type flags): {len(usable)}")
    for _, r in usable.iterrows():
        print(f"    {r['accession']:12} {str(r['n_samples']):>4} samples  "
              f"{r['inferred_cancer']:10} {r['title'][:60]}")
    if len(usable):
        print(f"\n  Next: python scripts/download_geo.py download "
              f"{usable.iloc[0]['accession']}")
        print(f"  Or let the pipeline do everything: "
              f"python scripts/download_geo.py {short}")


# ── Download ─────────────────────────────────────────────────────────────

def download_suppl(gse_id, raw_dir, files):
    """Fetch each supplementary file by its own URL.

    Not a recursive crawl: ftp.ncbi.nlm.nih.gov serves robots.txt with
    "Disallow: /", so wget -r refuses to follow the index links and comes
    back with nothing. The file list is already known from the inspection
    step, so each file is requested directly, which is also deterministic
    and lets a failure be reported per file.
    """
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/suppl/"
    failed = []
    for f in files:
        name = f["file"]
        out = os.path.join(raw_dir, name)
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"      {name} exists ({human_size(os.path.getsize(out))}) — skip")
            continue
        print(f"      {name} ({f['size']})")
        rc = subprocess.run(["wget", "-c", "-q", "--show-progress",
                             "-O", out, base + name]).returncode
        if rc != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
            failed.append(name)
            if os.path.exists(out) and os.path.getsize(out) == 0:
                os.remove(out)
    return failed


def sample_metadata(gse_id, out_dir):
    """Per-sample annotations via GEOparse.

    For sequencing series the GSM records normally carry no expression
    table — the counts live in the supplementary files — so this collects
    the annotations only, which is what the tables are reliable for.
    """
    try:
        import GEOparse
    except ImportError:
        print("    GEOparse not installed — skipping sample metadata")
        return None

    tmp_dir = os.path.join(out_dir, "_tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=tmp_dir, silent=True)
        rows = []
        for gsm_id, gsm in gse.gsms.items():
            row = {"gsm_id": gsm_id}
            for key, vals in gsm.metadata.items():
                row[key] = "; ".join(map(str, vals)) if isinstance(vals, list) else str(vals)
            rows.append(row)
        df = pd.DataFrame(rows)
        path = os.path.join(out_dir, "sample_metadata.tsv")
        df.to_csv(path, sep="\t", index=False)
        print(f"    Sample metadata: {len(df)} samples x {len(df.columns)} fields")

        with_table = sum(1 for g in gse.gsms.values()
                         if g.table is not None and not g.table.empty)
        print(f"    GSMs carrying an expression table: {with_table}/{len(gse.gsms)}"
              f"{' (normal for sequencing data)' if with_table == 0 else ''}")
        return df
    except Exception as e:
        print(f"    ERROR reading series: {e}")
        return None
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def download_gse(gse_id, out_root=None, quiet_inspect=False, ambiguity=""):
    if not quiet_inspect:
        print(f"\n{'=' * 60}")
        print(f"  Downloading {gse_id}")
        print(f"{'=' * 60}")
        print("\n  [1/4] Inspecting before download...")
    row, files = inspect_gse(gse_id, quiet=quiet_inspect)
    if row["has_raw"] == "no":
        print("\n  WARNING: no raw counts detected. Downloading anyway, but this"
              "\n  series may not be usable for differential expression.")

    # Default target when called on its own: the inferred cancer's folder
    out_dir = os.path.join(out_root or cancer_dir(row["inferred_cancer"]), gse_id)
    raw_dir = os.path.join(out_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    print("\n  [2/4] Sample annotations...")
    sample_metadata(gse_id, out_dir)

    print(f"\n  [3/4] Supplementary files ({len(files)})...")
    if files:
        pd.DataFrame(files).to_csv(os.path.join(out_dir, "suppl_files.csv"), index=False)
        # A <GSE>_RAW.tar only repackages the per-sample files, so when the
        # series also publishes a count matrix the bundle is redundant and
        # usually enormous: GSE162739 ships 7.5 GB of bundle next to a 1 MB
        # Raw_counts.tsv.gz holding what we actually need.
        wanted = files
        if any(f["kind"] == "raw" for f in files):
            bundles = [f for f in files if f["kind"] == "bundle"]
            if bundles:
                saved = sum(parse_size(f["size"]) for f in bundles)
                print(f"    Skipping {len(bundles)} redundant bundle(s), "
                      f"{human_size(saved)}: a count matrix is published separately")
                wanted = [f for f in files if f["kind"] != "bundle"]
        failed = download_suppl(gse_id, raw_dir, wanted)
    else:
        failed = []
        print("    none listed — the data may only be in SRA")
    got = sorted(f for f in os.listdir(raw_dir) if not f.startswith("."))
    on_disk = sum(os.path.getsize(os.path.join(raw_dir, f)) for f in got)
    print(f"    Files on disk: {len(got)} ({human_size(on_disk)})")
    if failed:
        print(f"    FAILED: {', '.join(failed)}")
        print(f"    Download by hand from {row['url']}")

    print("\n  [4/4] Writing summary...")
    rel = os.path.relpath(out_dir, os.path.dirname(DATASET_DIR))
    with open(os.path.join(out_dir, "summary.md"), "w") as f:
        f.write(f"# {gse_id} — {row['title']}\n\n")
        if ambiguity:
            # The filters could not settle this one, so say so where the
            # data is, with the command to remove it.
            f.write("## NEEDS REVIEW — this dataset is ambiguous\n\n")
            f.write("The automatic filters could not decide whether this series "
                    "belongs in the\nanalysis. It was downloaded so you can look "
                    "at it, not because it passed.\n\n")
            f.write(f"**Why:** {ambiguity}\n\n")
            f.write(f"Check the samples on the GEO page below. If this dataset is "
                    f"not needed,\ndelete it:\n\n```bash\nrm -rf {rel}\n```\n\n")
            f.write("---\n\n")
        f.write(f"- **Inferred cancer:** {row['inferred_cancer']}\n")
        f.write(f"- **Samples:** {row['n_samples']}\n")
        f.write(f"- **Platform:** {row['platform']}\n")
        f.write(f"- **Flags:** {row['flags']}\n")
        f.write(f"- **Raw counts available:** {row['has_raw']}\n")
        f.write(f"- **Verdict:** {row['verdict']}\n")
        f.write(f"- **Review status:** {'AMBIGUOUS — ' + ambiguity if ambiguity else 'passed every filter'}\n\n")
        f.write(f"## Supplementary files\n\n")
        for fl in files:
            f.write(f"- `{fl['file']}` ({fl['size']}) — {fl['kind']}\n")
        f.write(f"\n## Source\n{row['url']}\n")
    print(f"    summary.md written{' (marked AMBIGUOUS)' if ambiguity else ''}")
    print(f"\n  Done: {out_dir}")


def run_pipeline(target, max_download=None, include_flagged=False,
                 min_samples=10):
    """Search, inspect and download in one run, as TCGA and GTEx do.

    Only datasets that clear both filters are downloaded: raw counts present
    AND no sample-type flag. The two are independent — a series can ship raw
    counts and still be cell lines — so passing one is not enough. Flagged
    series are still searched, inspected and recorded; they are just not
    fetched, because the flags come from free-text metadata and deserve a
    human look before spending disk on them.
    """
    cancers = (CANCER_TYPES if target == "ALL" else
               {k: v for k, v in CANCER_TYPES.items() if v["short"] == target})
    if not cancers:
        valid = ", ".join(v["short"] for v in CANCER_TYPES.values())
        print(f"Unknown cancer: {target}. Valid: {valid}, all")
        sys.exit(1)

    label = "all 6 GI cancers" if target == "ALL" else target
    print(f"\n{'=' * 60}")
    print(f"  CLAUD-IA WP1 — GEO pipeline: {label}")
    print(f"{'=' * 60}")

    # ── [1/5] Search ────────────────────────────────────────────────
    print(f"\n  [1/5] Searching GEO...")
    found = []
    for project_id, info in cancers.items():
        print(f"    {info['label']}")
        results = search_geo_rnaseq(info["keywords"])
        for r in results:
            r["search_cancer"] = info["short"]
        found.extend(results)
        time.sleep(1)

    if not found:
        print("    No results.")
        return
    search_df = pd.DataFrame(found).drop_duplicates(subset="accession")
    for short, part in search_df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_dir(short), "search_results.csv"), index=False)
    print(f"    {len(search_df)} unique candidates")
    for flag, n in search_df["flags"].value_counts().items():
        print(f"      {n:4d}  {flag}")

    # ── [2/5] Inspect ───────────────────────────────────────────────
    print(f"\n  [2/5] Inspecting candidates (checking for raw counts)...")
    insp_csv = os.path.join(GEO_ROOT, "inspected.csv")
    already = set()
    if os.path.exists(insp_csv):
        already = set(pd.read_csv(insp_csv)["accession"])
    todo = [a for a in search_df["accession"] if a not in already]
    print(f"    {len(already & set(search_df['accession']))} already inspected, "
          f"{len(todo)} to go")
    for i, acc in enumerate(todo, start=1):
        time.sleep(NCBI_PAUSE)
        try:
            row, _ = inspect_gse(acc, quiet=True)
            print(f"    [{i}/{len(todo)}] {acc:12} raw={row['has_raw']:3} "
                  f"{row['flags'][:22]:22} {row['verdict']}")
        except Exception as e:
            print(f"    [{i}/{len(todo)}] {acc}: ERROR {e}")

    insp = pd.read_csv(insp_csv) if os.path.exists(insp_csv) else pd.DataFrame()
    insp = insp[insp["accession"].isin(set(search_df["accession"]))]

    # ── [3/5] Select, with a reason recorded for every candidate ────
    print(f"\n  [3/5] Deciding what to download...")
    if insp.empty:
        print("    Nothing inspected — stopping")
        return

    # Cheap filters first, from what the inspection already knows. Only the
    # candidates that survive get their files opened, because that costs
    # extra requests.
    # Two kinds of reason. A hard one is a fact: no raw counts, counts that
    # are not whole, a declared cell line. A soft one is a suspicion that
    # only a human can settle, so those are downloaded and marked ambiguous
    # in their summary.md rather than discarded.
    decisions = []
    to_verify = []
    for _, r in insp.iterrows():
        reasons, soft = [], []
        if r["has_raw"] != "yes":
            reasons.append(f"no raw counts: {str(r['verdict']).split(' — ')[0]}")
        flags = str(r.get("flags", ""))
        if flags and flags != "OK" and not include_flagged:
            reasons.append(f"sample-type flag: {flags}")
        # Sample count as a proxy for study design. A 6-sample series with a
        # 3-vs-3 layout is a mechanistic experiment, and those run on cell
        # lines far more often than on biopsies — but the abstract rarely
        # says so, so no text filter catches them. Next to TCGA's 1817
        # samples a 6-sample series adds little, while quietly importing
        # cell-line data costs a lot. --min-samples 0 turns this off.
        try:
            n = int(r.get("n_samples") or 0)
        except (TypeError, ValueError):
            n = 0
        if min_samples and 0 < n < min_samples:
            soft.append(f"only {n} samples (< {min_samples}): small studies are "
                        f"usually mechanistic experiments on cell lines")
        inferred = str(r.get("inferred_cancer", ""))
        if inferred in ("UNKNOWN", "AMBIGUOUS"):
            soft.append("the abstract names no cancer type clearly")
        elif inferred and inferred != target and target != "ALL":
            soft.append(f"text points at {inferred}, not {target}")
        row = {
            "accession": r["accession"],
            "cancer": r.get("inferred_cancer", ""),
            "n_samples": r.get("n_samples", ""),
            "size": human_size(r.get("download_bytes")
                               if pd.notna(r.get("download_bytes"))
                               else (r.get("suppl_bytes") or 0)),
            "reasons": reasons,
            "soft": soft,
            "title": str(r.get("title", ""))[:60],
        }
        (decisions if reasons else to_verify).append(row)

    # Content check: names and abstracts can both be wrong in the same
    # direction, towards a dataset looking more usable than it is.
    if to_verify:
        print(f"    Opening the files of {len(to_verify)} candidate(s) that"
              f" passed the name and abstract checks...")
    for row in to_verify:
        acc = row["accession"]
        files, _url = list_suppl_files(acc)
        chk = verify_contents(acc, files)
        if chk["raw_is_integer"] == "no":
            row["reasons"].append("content: 'raw' counts contain decimals, "
                                  "so they are not raw counts")
        if chk["content_flags"] and not include_flagged:
            row["reasons"].append(f"content: files mention {chk['content_flags']}")
        decisions.append(row)
        rec = dict(insp[insp["accession"] == acc].iloc[0])
        rec.update(chk)
        if row["reasons"]:
            rec["verdict"] = "REJECTED on content"
        record_inspection(rec, quiet=True)

    def verdict_of(d):
        if d["reasons"]:
            return "SKIP", "; ".join(d["reasons"])
        if d["soft"]:
            return "AMBIGUOUS", "; ".join(d["soft"])
        return "DOWNLOAD", "passes every filter"

    rows = []
    for d in decisions:
        decision, reason = verdict_of(d)
        rows.append({**d, "decision": decision, "reason": reason,
                     "ambiguity": "; ".join(d["soft"])})
    dec_df = pd.DataFrame(rows).drop(columns=["reasons", "soft"])
    dec_df = dec_df.sort_values(["decision", "accession"])
    cancer_of = dict(zip(search_df["accession"], search_df["search_cancer"]))
    dec_df["search_cancer"] = dec_df["accession"].map(cancer_of)
    for short, part in dec_df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_dir(short), "decisions.csv"), index=False)

    marks = {"DOWNLOAD": "KEEP", "AMBIGUOUS": " ?  ", "SKIP": "skip"}
    print(f"\n    Decision for each of the {len(dec_df)} candidates:")
    for _, r in dec_df.iterrows():
        print(f"      [{marks[r['decision']]}] {r['accession']:12} "
              f"{r['size']:>9}  {r['reason']}")

    print(f"\n    Why candidates were skipped:")
    skipped = dec_df[dec_df["decision"] == "SKIP"]
    for reason, n in skipped["reason"].str.split(";").str[0].value_counts().items():
        print(f"      {n:4d}  {reason.strip()}")

    amb = dec_df[dec_df["decision"] == "AMBIGUOUS"]
    if len(amb):
        print(f"\n    {len(amb)} marked AMBIGUOUS — downloaded, but review them:")
        for _, r in amb.iterrows():
            print(f"      {r['accession']:12} {r['reason']}")

    selected = dec_df[dec_df["decision"].isin(["DOWNLOAD", "AMBIGUOUS"])]
    if max_download and len(selected) > max_download:
        selected = selected.head(max_download)
        print(f"    capped at {max_download} by --max-download")

    if selected.empty:
        print(f"\n    Nothing passed every filter — nothing to download.")
        print(f"    Reasons are in dataset/GEO/<CANCER>/decisions.csv")
        print(f"    Options: review the UNCLEAR/BUNDLE ones by hand, or use"
              f" --include-flagged")
        print(f"    to fetch the flagged series anyway.")
        return

    total = sum(parse_size(x.replace(" ", "")) for x in selected["size"])
    print(f"\n    Downloading {len(selected)} dataset(s), {human_size(total)}:")
    for i, (_, r) in enumerate(selected.iterrows(), start=1):
        print(f"      {i:3}. {r['accession']:12} {str(r['n_samples']):>5} samples  "
              f"{r['size']:>9}  {r['cancer']:9} {r['title'][:46]}")


    # ── [4/5] Download ──────────────────────────────────────────────
    print(f"\n  [4/5] Downloading {len(selected)} datasets...")
    amb_of = dict(zip(dec_df["accession"], dec_df["ambiguity"]))
    for i, acc in enumerate(selected["accession"], start=1):
        tag = "  [AMBIGUOUS]" if amb_of.get(acc) else ""
        print(f"\n    --- [{i}/{len(selected)}] {acc}{tag} ---")
        try:
            download_gse(acc, out_root=cancer_dir(cancer_of.get(acc, target)),
                         quiet_inspect=True, ambiguity=amb_of.get(acc) or "")
        except Exception as e:
            print(f"    ERROR downloading {acc}: {e}")

    # ── [5/5] Summary ───────────────────────────────────────────────
    print(f"\n  [5/5] Summary")
    downloaded = [d for root, dirs, _ in os.walk(GEO_ROOT) for d in dirs
                  if d.startswith("GSE")]
    print(f"    Datasets on disk: {len(downloaded)}")
    n_raw = int((insp["has_raw"] == "yes").sum())
    n_clean = int((dec_df["decision"] == "DOWNLOAD").sum())
    print(f"    Funnel: {len(search_df)} found -> {len(insp)} inspected -> "
          f"{n_raw} with raw counts -> {n_clean} clean + {len(amb)} ambiguous")
    print(f"    Per-candidate reasons: dataset/GEO/<CANCER>/decisions.csv")
    print(f"\n  Next: python scripts/summary_report.py --md5")


# ── Main ─────────────────────────────────────────────────────────────────

CANCER_SHORTS = [v["short"] for v in CANCER_TYPES.values()]


def usage():
    shorts = "|".join(CANCER_SHORTS)
    print("Usage:")
    print(f"  python scripts/download_geo.py <{shorts}|all> [options]")
    print("        full pipeline: search -> inspect -> select -> download -> summary")
    print()
    print("  Options:")
    print("    --max-download N   cap how many datasets get downloaded")
    print("    --include-flagged  also download cell-line/organoid/single-cell series")
    print("    --min-samples N    skip series with fewer than N samples (default 10,")
    print("                       0 disables): small studies are usually cell lines")
    print()
    print("  Individual steps, for granular work:")
    print(f"    search <{shorts}|all>")
    print("    inspect <GSE_ID>")
    print(f"    inspect all <{shorts}> [--ok-only] [--recheck]")
    print("    download <GSE_ID>")
    sys.exit(1)


def run_search(target):
    if target == "ALL":
        cancers = CANCER_TYPES
    else:
        cancers = {k: v for k, v in CANCER_TYPES.items() if v["short"] == target}
        if not cancers:
            valid = ", ".join(v["short"] for v in CANCER_TYPES.values())
            print(f"Unknown cancer: {target}. Valid: {valid}")
            sys.exit(1)

    all_results = []
    for project_id, info in cancers.items():
        print(f"\n  Searching GEO: {info['label']}")
        found = search_geo_rnaseq(info["keywords"])
        for r in found:
            r["search_cancer"] = info["short"]
        all_results.extend(found)
        time.sleep(1)

    if not all_results:
        print("\n  No results.")
        return

    df = pd.DataFrame(all_results).drop_duplicates(subset="accession")
    for short, part in df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_dir(short), "search_results.csv"), index=False)
    out_path = os.path.join(GEO_ROOT, "<CANCER>", "search_results.csv")

    print(f"\n{'=' * 60}")
    print(f"  {len(df)} unique datasets")
    print(f"{'=' * 60}\n")
    print(df[["accession", "search_cancer", "inferred_cancer",
              "n_samples", "flags"]].to_string(index=False))
    print(f"\n  Flag breakdown:")
    for flag, n in df["flags"].value_counts().items():
        print(f"    {n:4d}  {flag}")
    print(f"\n  Saved: dataset/GEO/<CANCER>/search_results.csv")
    n_ok = int((df["flags"] == "OK").sum())
    print(f"\n  Next: inspect them in one go (checks for raw counts):")
    for short in sorted(df["search_cancer"].unique()):
        print(f"    python scripts/download_geo.py inspect all {short}")
    print(f"    ({n_ok} of {len(df)} carry no sample-type flag; add --ok-only"
          f" to limit to those)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()

    mode = sys.argv[1].lower()
    # The pipeline form takes a cancer code and no sub-command
    if mode.upper() not in CANCER_SHORTS and mode.upper() != "ALL" and len(sys.argv) < 3:
        usage()

    start_logging()
    arg = sys.argv[2] if len(sys.argv) > 2 else ""

    if mode == "search":
        run_search(arg.upper())
    elif mode.upper() in CANCER_SHORTS or mode.upper() == "ALL":
        # No sub-command: the whole pipeline for that cancer
        run_pipeline(mode.upper(),
                     max_download=(int(sys.argv[sys.argv.index("--max-download") + 1])
                                   if "--max-download" in sys.argv else None),
                     include_flagged="--include-flagged" in sys.argv,
                     min_samples=(int(sys.argv[sys.argv.index("--min-samples") + 1])
                                  if "--min-samples" in sys.argv else 10))
    elif mode == "inspect":
        if arg.lower() == "all":
            if len(sys.argv) < 4 or sys.argv[3].upper() not in CANCER_SHORTS:
                print(f"  Which cancer? e.g. inspect all CHOL")
                sys.exit(1)
            inspect_batch(sys.argv[3].upper(),
                          ok_only="--ok-only" in sys.argv,
                          recheck="--recheck" in sys.argv)
        else:
            inspect_gse(arg.upper())
    elif mode == "download":
        download_gse(arg.upper())
    else:
        usage()
