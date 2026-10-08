"""
download_geo.py — CLAUD-IA WP1: find and download bulk RNA-seq from GEO.

Usage:
    python scripts/download_geo.py CHOL          # whole pipeline for one cancer
    python scripts/download_geo.py all           # the six GI cancers
    python scripts/download_geo.py inspect GSE107943    # one series, in detail
    python scripts/download_geo.py search CHOL          # candidates only

Options:
    --max-download N    stop after N datasets
    --include-flagged   fetch the cell-line/organoid series too
    --min-samples N     suspicion threshold (default 10, 0 disables)

Why this is curation and not a download
---------------------------------------
TCGA and GTEx each publish one uniform matrix: copy it and you are done. GEO
is thousands of laboratories, each uploading what they liked, and most of
what they upload is not what WP1 needs. Of 204 cholangiocarcinoma candidates,
5 held bulk RNA-seq counts of human tissue. So the work is to reject, with a
stated reason for every rejection.

What "clean" means here
-----------------------
Clean requires positive evidence, not the absence of detected problems. That
distinction is the whole design. The earlier version passed anything no
filter had flagged, which quietly admitted series whose authors never said
what their samples were, and series whose files could not be opened to check.
GSE254942 is the cautionary case: 106 samples, Ensembl ids, whole-number
counts, every filter satisfied — and its sample records say "cell line: 007b".

So a series is only clean when the submitter's own per-sample fields declare
tissue. Anything unconfirmed is marked ambiguous, downloaded, and flagged for
review in its own summary.md. Hard facts reject; suspicions mark.

The four gates, cheapest first
------------------------------
    1. text      the abstract and title              1 request
    2. names     the supplementary file names        1 request
    3. contents  the first KB of the count matrix    1-3 requests
    4. samples   what the submitter declared         1 request

Gates 1 and 2 run over every candidate and are cached in inspected.csv.
Gates 3 and 4 only run on the survivors, because they cost more; gate 4 is
the decisive one and reads rather than infers.

Output:
    dataset/GEO/
    ├── inspected.csv               gates 1-2 per candidate, cached
    └── <CANCER>/
        ├── search_results.csv      what the search found
        ├── decisions.csv           verdict and reason for every candidate
        └── <GSE>/
            ├── summary.md          the dataset's card, NEEDS REVIEW if ambiguous
            ├── sample_metadata.tsv per-sample annotations
            ├── suppl_files.csv     files and their classification
            └── raw/                the supplementary files
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
from config import (DATASET_DIR, CURATION_DIR, CANCER_TYPES,
                    CANCER_ORGAN_PATTERNS)
from run_logger import start_logging

# ══════════════════════════════════════════════════════════════════════════
#  Configuration
# ══════════════════════════════════════════════════════════════════════════

GEO_ROOT = os.path.join(DATASET_DIR, "GEO")       # the files
GEO_REC = os.path.join(CURATION_DIR, "GEO")       # why we kept them
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACC_CGI = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# NCBI allows three requests a second without an API key
NCBI_PAUSE = 0.35

CANCER_SHORTS = [v["short"] for v in CANCER_TYPES.values()]

# ── Sample-type patterns ─────────────────────────────────────────────────
# Matched against abstracts, file names and column headers. Every entry was
# added because a real series slipped past without it; the comments say which.
FLAG_PATTERNS = {
    "CELL_LINE": (
        r"cell\s*lines?|ccle|"
        # named lines that recur in GI work
        r"\bhek-?293\b|\bhela\b|panc-?1|mia\s*paca|bxpc-?3|capan|aspc-?1|"
        r"\bhct-?116\b|\bhepg2\b|\bsw480\b|caco-?2|\ba549\b|mcf-?7|\bu-?87\b|"
        r"\bu-?138\b|\blo2\b|\bthle-?2?\b|hucct|\brbe\b|qbc939|huh-?\d+|"
        r"tfk-?1|\bccl?p-?1\b|snu-?\d{3,4}|egi-?1|mz-?cha|gbc-?sd|\bnoz\b|"
        r"kyse|\bte-?\d{1,2}\b|\beca?-?109\b|\boe-?(19|21|33)\b|flo-?1|"
        r"kato|\bags\b|mkn-?\d+|sgc-?7901|bgc-?823|hgc-?27|nci-?n87|"
        # engineered cultures. A digit is required so that "short" and
        # "should" do not match: GSE261202 ships shIGF2BP3.gene_count.csv.gz
        r"\bshrna\b|\bsirna\b|knock-?down|knock-?out|\bko\b|\bkd\b|\bcrispr\b|"
        r"transfect|overexpress|\bparental\b|\bclones?\b|[a-z0-9]{2,}ko\b|"
        r"\bsh\d\b|\bsi\d\b|\bsh[a-z]{2,}\d[a-z0-9]*\b|\bscr\d?\b|scramble|"
        r"\bshctrl\b|\bshnc\b|\bshluc\b|\bnegative control\b"
    ),
    "ORGANOID": r"organoid|\bspheroid\b",
    "XENOGRAFT": r"xenograft|\bpdx\b",
    # HFD (high-fat diet) studies are mouse work
    "NON_HUMAN": (r"\bmouse\b|\bmice\b|murine|mus musculus|\brat\b|rattus|"
                  r"zebrafish|drosophila|\bhfd\b|high[-\s]?fat diet"),
    # Not bulk. The file names give these away even when the abstract does not
    "SINGLE_CELL": (r"single[-\s]?cell|single[-\s]?nucle|\bscrna\b|\bsnrna\b|"
                    r"10x genomics|fragments\.tsv|barcodes\.tsv|features\.tsv|"
                    r"matrix\.mtx|cell[-_\s]?annotation"),
    # Not RNA-seq at all
    "NOT_RNASEQ": (r"\batac\b|\bchip[-\s]?seq\b|\bcut&?run\b|bisulfite|"
                   r"\bwgbs\b|\bmethylat"),
    # Liquid biopsy is not a tissue biopsy: GSE183635 is 2351 platelet samples
    "NOT_TISSUE": (r"platelets?|liquid biops|\bplasma\b|\bserum\b|whole blood|"
                   r"\bpbmc\b|buffy coat|cell[-\s]?free|\bcfrna\b|\bctdna\b|"
                   r"exosom|extracellular vesicle|\bsaliva\b|\burine\b"),
}

# ── Supplementary file classification ────────────────────────────────────
# Order matters: a file called "Normcount" is normalised despite the word
# "count", and "circRNAs_count" is not gene level despite both words.
NON_GENE_MARKERS = (r"\bcircrnas?\b|circular rna|\bisoforms?\b|\btranscripts?\b|"
                    r"\bmirnas?\b|microrna|\bexons?\b|\bjunctions?\b|\bpeaks?\b|"
                    r"\bpromoters?\b|\bsplic|\blncrnas?\b|long non-?coding|"
                    r"\bpirnas?\b|\bsnornas?\b|\btrnas?\b")
NORMALISED_MARKERS = (r"fpkm|rpkm|\btpm\b|\bcpm\b|norm|scaled|vst|rlog|"
                      r"z[-_]?score|deseq|edger|\blog2\b")
RAW_MARKERS = r"raw|count|htseq|featurecount|star|\bexpected_count\b|quant"
# "<GSE>_RAW.tar" is GEO's generic container for per-sample files, so RAW
# there means "the raw files", not "raw counts": one held scATAC fragments
BUNDLE_MARKERS = r"_raw\.(tar|zip)$"

# ── What submitters write in the per-sample fields ───────────────────────
TISSUE_WORDS = (r"\btissues?\b|\bbiops|\bresect|surgical|\bpatients?\b|\bffpe\b|"
                r"\bspecimens?\b|adjacent normal|\bparaffin\b")
CULTURE_WORDS = (r"cell\s*lines?|\bcells?\b|organoid|spheroid|xenograft|\bpdx\b|"
                 r"passage|cultur|\bmedium\b|\bdmem\b|\brpmi\b|\bfbs\b")
# "squamous cell carcinoma" contains the word "cell", which counted as
# evidence of culture and tipped GSE248518 — oesophageal tumour tissues —
# into the culture bucket. Diagnoses are stripped of it before matching.
PATHOLOGY_CELL = [
    (r"\b(squamous|clear|small|large|non-small|transitional|germ|spindle)\s+cell\b", r"\1"),
    (r"\bcell\s+(carcinomas?|tumou?rs?|cancers?|neoplasms?)\b", r"\1"),
]


def strip_pathology(text):
    for pat, rep in PATHOLOGY_CELL:
        text = re.sub(pat, rep, text)
    return text


def looks_like_cell_line(value):
    """Is this the name of a cell line, or a sentence in the wrong field?

    GSE248518 records treatment response under "cell line:", writing
    "sensitivity to chemotherapy" where a line name belongs. A name is one
    or two tokens and usually carries a digit — 007b, flo1, oe33, hucct1 —
    whereas a misused field holds an English phrase.
    """
    if flag_sample_issues(value):
        return True
    words = value.split()
    return len(words) <= 2 and bool(re.search(r"\d", value))
# Organs outside WP1: a search for oesophageal cancer returned GSE253064,
# whose samples declare "human urinary bladder tumor"
OTHER_ORGANS = (r"\bbladder\b|\bbreast\b|\blung\b|\bprostate\b|\bglioma\b|"
                r"glioblastoma|\bbrain\b|\bkidney\b|\brenal\b|\bovar|\bcervi|"
                r"\bthyroid\b|melanoma|leukemi|leukaemi|lymphoma|\bmyeloma\b|"
                r"head and neck|\bskin\b|\buterus\b|endometri|\btestic")

# ── Gene identifiers, to tell a count matrix from something else ─────────
GENE_ID_PATTERNS = {
    "ensembl-human": r"^ENSG\d{11}",
    "ensembl-mouse": r"^ENSMUSG\d{11}",
    "ensembl-transcript": r"^ENS[A-Z]*T\d{8}",
    "refseq": r"^(NM_|NR_|XM_)\d+",
    "symbol": r"^[A-Z][A-Z0-9%s-]{1,20}$" % "_",
}


def cancer_dir(short):
    """dataset/GEO/<CANCER>/ — where the downloaded files go."""
    path = os.path.join(GEO_ROOT, short.upper())
    os.makedirs(path, exist_ok=True)
    return path


def cancer_rec(short):
    """curation/GEO/<CANCER>/ — where the record of the decisions goes.

    Separate from the data because dataset/ is gitignored wholesale: the
    files can be downloaded again, the reasoning cannot. The verdict cache
    sits at curation/GEO/, not under a cancer, because a series' verdict
    does not depend on which search found it."""
    path = os.path.join(GEO_REC, short.upper())
    os.makedirs(path, exist_ok=True)
    return path


def tokenise(text):
    """Lowercase and split on punctuation.

    Needed everywhere a pattern meets a file name: regex counts "_" as a word
    character, so \\blo2\\b never matches GSE276342_LO2_raw_count.txt.gz.
    """
    return re.sub(r"[^a-z0-9&]+", " ", text.lower())


def parse_size(text):
    """'1.8M' -> bytes, as the HTTPS directory listing writes sizes."""
    m = re.match(r"([0-9.]+)\s*([KMGT]?)", (text or "").strip(), re.I)
    if not m:
        return 0
    mult = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}
    return int(float(m.group(1)) * mult[m.group(2).upper()])


def human_size(n):
    n = float(n or 0)
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024


# ══════════════════════════════════════════════════════════════════════════
#  NCBI access
# ══════════════════════════════════════════════════════════════════════════

def esummary(gds_id):
    r = requests.get(f"{EUTILS}/esummary.fcgi",
                     params={"db": "gds", "id": gds_id, "retmode": "json"},
                     timeout=60)
    r.raise_for_status()
    return r.json()["result"].get(gds_id, {})


def fetch_series(gse_id):
    """The series record for one accession."""
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gds", "term": f"{gse_id}[ACCN]",
                             "retmode": "json"}, timeout=60)
    r.raise_for_status()
    ids = r.json()["esearchresult"].get("idlist", [])
    if not ids:
        raise ValueError(f"No NCBI record for {gse_id}")
    time.sleep(NCBI_PAUSE)
    return esummary(ids[0])


def search_candidates(keywords):
    """Every bulk RNA-seq series for one cancer. No cap, on purpose.

    NCBI returns the most recent first, so a capped search is a biased
    sample rather than a smaller one: what gets published now is mostly
    single-cell and cell-line work, while bulk RNA-seq of biopsies is
    commoner among the older entries.
    """
    query = (f'({" OR ".join(keywords)}) AND '
             '"expression profiling by high throughput sequencing"[DataSet Type] AND '
             '"homo sapiens"[Organism]')
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gds", "term": query, "retmax": 100000,
                             "retmode": "json"}, timeout=120)
    r.raise_for_status()
    res = r.json()["esearchresult"]
    ids = res.get("idlist", [])
    print(f"    {res.get('count', 0)} datasets match; retrieving all {len(ids)}")
    mins = len(ids) * 2 * NCBI_PAUSE / 60
    if mins > 2:
        print(f"    ~{mins:.0f} min to list and inspect "
              f"(NCBI allows 3 requests/second). Resumable.")

    out = []
    for gds_id in ids:
        time.sleep(NCBI_PAUSE)
        meta = esummary(gds_id)
        if not str(meta.get("accession", "")).startswith("GSE"):
            continue
        out.append(meta)
    return out


def list_suppl_files(gse_id):
    """Names, sizes and classification of a series' supplementary files.

    NCBI serves its FTP tree over HTTPS too, which is faster and less
    likely to be blocked.
    """
    url = f"{FTP_HTTPS}/{gse_id[:-3]}nnn/{gse_id}/suppl/"
    try:
        r = requests.get(url, timeout=90)
    except Exception:
        return [], url
    if not r.ok:
        return [], url
    files = []
    for m in re.finditer(
            r'href="([^"?/][^"]*)"[^<]*</a>\s*([0-9-]{10}\s[0-9:]{5})?\s*([0-9.]+[KMG]?)?',
            r.text):
        name = m.group(1)
        if name.startswith("http") or name == "../":
            continue
        files.append({"file": name, "size": (m.group(3) or "").strip(),
                      "kind": classify_suppl(name)})
    return files, url


def peek_file(gse_id, fname, max_bytes=262144):
    """The first few KB of a supplementary file, over an HTTP range request.

    Decompresses gzip on the fly. Returns "" for anything it cannot read as
    text, which includes .xlsx and .RData: a zip or an R binary cannot be
    interpreted from a fragment, and that failure has to count as unknown
    rather than as a pass.
    """
    url = f"{FTP_HTTPS}/{gse_id[:-3]}nnn/{gse_id}/suppl/{fname}"
    try:
        r = requests.get(url, headers={"Range": f"bytes=0-{max_bytes}"},
                         timeout=90)
        if r.status_code not in (200, 206):
            return ""
        data = r.content
        if fname.endswith(".gz"):
            import zlib
            data = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data)
        text = data.decode("utf-8", "replace")
        # A binary that happens to decode is still not a table
        if text.count("�") > len(text) * 0.02:
            return ""
        return text
    except Exception:
        return ""


def fetch_sample_records(gse_id):
    """Every GSM's source name and characteristics, in one plain-text request.

    These two fields are the only place in GEO where sample type is declared
    rather than implied, which makes this the decisive gate.
    """
    try:
        r = requests.get(ACC_CGI, params={"acc": gse_id, "targ": "gsm",
                                          "form": "text", "view": "brief"},
                         timeout=120)
        if not r.ok:
            return []
    except Exception:
        return []
    values = []
    for line in r.text.splitlines():
        line = line.strip()
        if line.startswith(("!Sample_source_name", "!Sample_characteristics")):
            values.append(line.split("=", 1)[1].strip().lower()
                          if "=" in line else line.lower())
    return values


# ══════════════════════════════════════════════════════════════════════════
#  Classification helpers
# ══════════════════════════════════════════════════════════════════════════

def classify_suppl(fname):
    """raw / normalised / non_gene / bundle / other, from the name alone."""
    low = fname.lower()
    if re.search(BUNDLE_MARKERS, low):
        return "bundle"
    tokens = tokenise(low)
    if re.search(NON_GENE_MARKERS, tokens):
        return "non_gene"
    if re.search(NORMALISED_MARKERS, tokens):
        return "normalised"
    if re.search(RAW_MARKERS, tokens):
        return "raw"
    return "other"


def infer_cancer(title, summary):
    """Which WP1 cancer a paper is about.

    Separate from the search keywords on purpose: those are phrases for
    NCBI, and matching them literally labelled "Esophageal Adenocarcinoma"
    as UNKNOWN because it does not contain the phrase "esophageal cancer".
    """
    text = f"{title} {summary}".lower()
    hits = [c for c, pat in CANCER_ORGAN_PATTERNS.items() if re.search(pat, text)]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        return "UNKNOWN"
    # Several organs named: the title outranks the abstract, and
    # CANCER_ORGAN_PATTERNS is ordered most specific first, so
    # "intrahepatic cholangiocarcinoma" resolves to CHOL and not LIHC.
    in_title = [c for c in hits if re.search(CANCER_ORGAN_PATTERNS[c], title.lower())]
    pool = set(in_title or hits)
    for c in CANCER_ORGAN_PATTERNS:
        if c in pool:
            return c
    return "AMBIGUOUS"


def flag_sample_issues(*texts):
    """Sample-type flags found in free text, file names or column headers."""
    raw = " ".join(t or "" for t in texts).lower()
    tokens = tokenise(raw)
    return [name for name, pat in FLAG_PATTERNS.items()
            if re.search(pat, raw) or re.search(pat, tokens)]


def split_fields(line):
    return [f.strip().strip('"').strip("'")
            for f in re.split(r"[\t,;]|\s{2,}|\s", line.strip()) if f != ""]


def values_are_integers(text):
    """Are the numbers whole? True / False / None if nothing was testable.

    Tokenising beats a regex over the raw text: fields can be separated by
    spaces as well as tabs, and the first field of each row is a gene id
    whose version suffix (ENSG00000240361.2) looks like a decimal.

    "5.0" counts as whole, because R writes integers that way.
    """
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    checked = 0
    for line in lines[1:60]:
        for v in split_fields(line)[1:]:
            if not re.fullmatch(r"-?\d+(\.\d+)?([eE][-+]?\d+)?", v):
                continue
            checked += 1
            if float(v) != int(float(v)):
                return False
    return True if checked else None


def gene_id_kind(text, max_cols=8):
    """Which identifier scheme the matrix uses, from its leading columns.

    Scans the first few columns rather than only the first, because the
    identifier is not always there: GSE107943 opens with a row number and
    puts Ensembl ids in column three, behind Chr.

    Also catches mouse data that called itself human, since the ids read
    ENSMUSG.
    """
    rows = [split_fields(l) for l in text.splitlines() if l.strip()][1:40]
    rows = [r for r in rows if r]
    if not rows:
        return ""
    for col in range(min(max_cols, max(len(r) for r in rows))):
        vals = [r[col] for r in rows if len(r) > col]
        if not vals:
            continue
        for kind, pat in GENE_ID_PATTERNS.items():
            if sum(bool(re.match(pat, v, re.I)) for v in vals) >= len(vals) * 0.7:
                return kind
    return "unrecognised"


def count_numeric_columns(text):
    """How many numeric columns the matrix has, i.e. how many samples.

    Compared against the sample count GEO declares. A mismatch means the
    file is not the matrix we assumed it was.
    """
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return 0
    best = 0
    for line in lines[1:10]:
        fields = split_fields(line)
        n = sum(1 for v in fields[1:]
                if re.fullmatch(r"-?\d+(\.\d+)?([eE][-+]?\d+)?", v))
        best = max(best, n)
    return best


# ══════════════════════════════════════════════════════════════════════════
#  The four gates
#
#  Each returns (hard, soft): hard reasons reject the series, soft reasons
#  mark it ambiguous. A gate never decides on its own; decide() combines them.
# ══════════════════════════════════════════════════════════════════════════

def gate_text(meta, target, min_samples):
    """Gate 1 — the title and abstract.

    Cheap and broad, but it only sees what the author chose to write, so it
    produces one hard reason (a declared cell line) and several suspicions.
    """
    title = str(meta.get("title", ""))
    summary = str(meta.get("summary", ""))
    hard, soft = [], []

    # A suspicion, not a fact. The abstract is the weakest evidence here and
    # it was taking the strongest action: plenty of tissue studies say they
    # validated in cell lines. Of 24 ESCA candidates rejected on the
    # abstract alone, GSE235537 declares "esophageal squamous carcinoma"
    # tissue across 22 samples — and it would have been the only confirmed
    # dataset for that cancer. Gate 4 reads what the samples are, so it
    # settles this one.
    flags = flag_sample_issues(title, summary)
    if flags:
        soft.append(f"the abstract mentions {', '.join(flags)}, which the "
                    f"sample records did not confirm")

    if re.search(r"circ(ular)?[-\s]?rna|\bmirna\b|microrna|\blncrna\b", title.lower()):
        soft.append("the title is about non-coding RNA, so gene-level counts "
                    "may be a by-product")

    cancer = infer_cancer(title, summary)
    if cancer in ("UNKNOWN", "AMBIGUOUS"):
        soft.append("the abstract names no cancer type clearly")
    elif target != "ALL" and cancer != target:
        soft.append(f"the text points at {cancer}, not {target}")

    try:
        n = int(meta.get("n_samples") or 0)
    except (TypeError, ValueError):
        n = 0
    if min_samples and 0 < n < min_samples:
        soft.append(f"only {n} samples (< {min_samples}): small studies are "
                    f"usually mechanistic experiments on cell lines")
    return hard, soft, {"inferred_cancer": cancer, "n_samples": n,
                        "abstract_flags": ", ".join(flags)}


def gate_files(gse_id, files):
    """Gate 2 — the supplementary file names.

    Often franker than the abstract: an author names files to work with
    them, not to publish them. GSE162739's abstract never says "cell line",
    but its files are Bap1KO-C91S-vs-PAR_diffexp.
    """
    hard, soft = [], []
    kinds = {f["kind"] for f in files}

    if not files:
        hard.append("no supplementary files: the data may only be in SRA")
    elif "raw" in kinds:
        pass
    elif "bundle" in kinds:
        hard.append("only a per-sample archive (<GSE>_RAW.tar), no count matrix")
    elif "non_gene" in kinds:
        hard.append("only circRNA/isoform/miRNA counts, not gene level")
    elif "normalised" in kinds:
        hard.append("only normalised values (FPKM/TPM/CPM), unusable by DESeq2")
    else:
        hard.append("no file recognisable as a count matrix")

    flags = flag_sample_issues(" ".join(f["file"] for f in files))
    if flags:
        hard.append(f"sample-type flag in the file names: {', '.join(flags)}")
    return hard, soft, {"has_raw": "yes" if "raw" in kinds else "no",
                        "n_suppl_files": len(files),
                        "file_flags": ", ".join(flags)}


def gate_contents(gse_id, files, declared_samples):
    """Gate 3 — the first KB of the count matrix.

    File names lie: GSE335154 ships Rawcount.csv.gz full of decimals, and
    GSE179443's space-separated Raw_gene_counts is not whole either. This
    opens the file and checks the values, the identifiers, the column count
    and the column names.
    """
    hard, soft = [], []
    raw = next((f for f in files if f["kind"] == "raw"), None)
    info = {"raw_is_integer": "", "gene_ids": "", "matrix_columns": ""}
    if not raw:
        return hard, soft, info

    text = peek_file(gse_id, raw["file"])
    if not text:
        ext = os.path.splitext(raw["file"])[1] or "an unknown format"
        soft.append(f"the counts are in {ext}, which cannot be read from a byte "
                    f"range: check by hand that the values are whole and the "
                    f"columns are tissue")
        return hard, soft, info

    whole = values_are_integers(text)
    info["raw_is_integer"] = {True: "yes", False: "no", None: ""}[whole]
    if whole is False:
        hard.append("the 'raw' counts contain decimals, so they are not counts")

    # The name of the identifier column says what the rows are. GSE248518
    # ships "all.count.txt.gz" whose first column is lncRNA_id, holding
    # NONCODE transcript ids: tissue samples, but not gene-level counts, and
    # no file name hinted at it.
    id_label = split_fields(text.splitlines()[0])[0] if text.splitlines() else ""
    if id_label and re.search(NON_GENE_MARKERS, tokenise(id_label)):
        hard.append(f"the identifier column is '{id_label}', so the rows are "
                    f"not genes")

    ids = gene_id_kind(text)
    info["gene_ids"] = ids
    if ids == "ensembl-mouse":
        hard.append("the gene identifiers are ENSMUSG: this is mouse data")
    elif ids in ("ensembl-transcript",):
        hard.append("the identifiers are transcripts, not genes")
    elif ids in ("", "unrecognised"):
        soft.append("the first column does not look like gene identifiers")

    ncol = count_numeric_columns(text)
    info["matrix_columns"] = ncol
    if declared_samples and ncol:
        if not (0.5 * declared_samples <= ncol <= 2 * declared_samples):
            soft.append(f"the matrix has {ncol} numeric columns but GEO declares "
                        f"{declared_samples} samples")

    # Column names describe the samples, and are the frankest metadata of
    # all: GSE312961's are GBC-SD-CASE and NOZ-CASE, two gallbladder lines
    header_flags = flag_sample_issues(text.splitlines()[0][:4000])
    if header_flags:
        hard.append(f"sample-type flag in the matrix header: "
                    f"{', '.join(header_flags)}")

    # Small annotation files, and a bundle's filelist, often name the lines
    for f in files:
        if f["kind"] == "bundle" or parse_size(f["size"]) > 2 * 1024 ** 2:
            continue
        if f["kind"] != "other":
            continue
        extra = flag_sample_issues(peek_file(gse_id, f["file"], 65536)[:20000])
        if extra:
            hard.append(f"sample-type flag inside {f['file']}: {', '.join(extra)}")
            break
    return hard, soft, info


def gate_samples(gse_id, target, include_flagged):
    """Gate 4 — what the submitter declared, per sample. The decisive one.

    Every gate above infers; this one reads. !Sample_source_name_ch1 and
    !Sample_characteristics_ch1 are the only fields in GEO where sample type
    is stated outright, and they settle cases nothing else could:

        GSE275730   cell line: cell line_TE-5        culture
        GSE254942   cell line: 007b                  culture
        GSE253064   tissue: human urinary bladder    wrong organ
        GSE107943   tissue: Tumor (intrahepatic CCA) tissue

    Clean requires reaching "tissue" here. Silence is not a pass.
    """
    hard, soft = [], []
    values = fetch_sample_records(gse_id)
    info = {"declared": "", "sample_detail": ""}
    if not values:
        soft.append("the sample records could not be read, so the sample type "
                    "is unconfirmed")
        return hard, soft, info

    # Split "label: value" so patterns only ever see the value. The label
    # itself must not be evidence: "cell line: sensitivity to chemotherapy"
    # is a misused field on tissue samples, not a cell line.
    pairs = []
    for v in values:
        label, _, val = v.partition(":")
        if val:
            pairs.append((label.strip(), strip_pathology(val.strip())))
        else:
            pairs.append(("", strip_pathology(v)))
    vals = [val for _, val in pairs if val]

    named_lines = [val for label, val in pairs
                   if re.match(r"cell\s*lines?$", label) and val
                   and val not in ("na", "n/a", "none", "not applicable", "-")
                   and looks_like_cell_line(val)]
    flags = set(flag_sample_issues(" ; ".join(vals)))
    if named_lines:
        flags.add("CELL_LINE")

    declared_tissue = [val for label, val in pairs
                       if re.match(r"(tissues?|organ|source|sample type)$", label)
                       and re.search(TISSUE_WORDS, val)]
    tissue_hits = sum(bool(re.search(TISSUE_WORDS, v)) for v in vals)
    culture_hits = sum(bool(re.search(CULTURE_WORDS, v)) for v in vals)
    organs = sorted({m.group(0) for v in vals
                     for m in [re.search(OTHER_ORGANS, v)] if m})

    # Order matters. A named line or an explicit flag settles it; otherwise
    # a "tissue:" field that names tissue outranks loose culture wording.
    if flags:
        info["declared"] = "culture"
    elif declared_tissue:
        info["declared"] = "tissue"
    elif culture_hits and culture_hits >= tissue_hits:
        info["declared"] = "culture"
    elif tissue_hits:
        info["declared"] = "tissue"
    elif (target in CANCER_ORGAN_PATTERNS
          and any(re.search(CANCER_ORGAN_PATTERNS[target], v) for v in vals)
          and not culture_hits):
        # Naming the organ with no hint of culture is tissue evidence too.
        # GSE291098 says "esophagus tumor", which carries none of the words
        # in TISSUE_WORDS but is plainly a biopsy.
        info["declared"] = "tissue"
    else:
        info["declared"] = "unstated"
    info["sample_detail"] = (named_lines or declared_tissue or vals or [""])[0][:120]

    # Another organ only disqualifies when our own is absent. A pancreatic
    # study with lung metastases names both, and metastatic samples are
    # wanted — TCGA code 06 is kept for the same reason. Without our organ
    # but with metastasis wording it could still be a metastasis of our
    # cancer, so that is a suspicion rather than a fact.
    own_organ = (target != "ALL" and target in CANCER_ORGAN_PATTERNS
                 and re.search(CANCER_ORGAN_PATTERNS[target], " ; ".join(vals)))
    if organs and not own_organ:
        if re.search(r"metasta", " ; ".join(vals)):
            soft.append(f"the samples name {', '.join(organs)} in a metastasis "
                        f"context but never {target}: check whether these are "
                        f"metastases of our cancer or another primary")
        else:
            hard.append(f"the samples declare another organ: {', '.join(organs)}")
    elif info["declared"] == "culture" and not include_flagged:
        detail = f" ({', '.join(sorted(flags))})" if flags else ""
        hard.append(f"the samples declare culture{detail}: {info['sample_detail']}")
    elif info["declared"] == "unstated":
        soft.append("the samples declare neither tissue nor cell line, so the "
                    "type could not be confirmed")
    return hard, soft, info


# ══════════════════════════════════════════════════════════════════════════
#  Inspection (gates 1-2) and decision (gates 3-4)
# ══════════════════════════════════════════════════════════════════════════

def inspect_gse(gse_id, target="ALL", min_samples=10, quiet=False):
    """Run the two cheap gates over one series and record the outcome.

    Cached in inspected.csv, because these two gates run over every
    candidate — 2080 of them for colorectal — and nothing in them depends
    on which cancer's search found the series.
    """
    meta = fetch_series(gse_id)
    files, url = list_suppl_files(gse_id)

    h1, s1, i1 = gate_text(meta, target, min_samples)
    h2, s2, i2 = gate_files(gse_id, files)

    row = {"accession": gse_id, **i1, **i2,
           "hard": " | ".join(h1 + h2), "soft": " | ".join(s1 + s2),
           "suppl_bytes": sum(parse_size(f["size"]) for f in files),
           # What a download would fetch: a bundle is skipped when a count
           # matrix is published next to it
           "download_bytes": sum(
               parse_size(f["size"]) for f in files
               if not (f["kind"] == "bundle"
                       and any(g["kind"] == "raw" for g in files))),
           "title": str(meta.get("title", ""))[:140],
           "platform": f"GPL{meta.get('gpl', '?')}",
           "url": f"{ACC_CGI}?acc={gse_id}"}

    if not quiet:
        print(f"\n{'=' * 60}\n  {gse_id}\n{'=' * 60}")
        print(f"\n  Title:    {row['title'][:90]}")
        print(f"  Samples:  {row['n_samples']}   Platform: {row['platform']}")
        print(f"  Cancer:   {row['inferred_cancer']}")
        print(f"\n  Supplementary files ({len(files)}):")
        for f in files or []:
            print(f"    [{f['kind']:10}] {f['size']:>8}  {f['file']}")
        if not files:
            print(f"    none listed — check {url}")
        for label, items in (("Rejected by", h1 + h2), ("Unclear", s1 + s2)):
            for x in items:
                print(f"  {label}: {x}")
        if not (h1 + h2):
            print("  Passes the text and file-name gates")
    record_inspection(row, quiet=quiet)
    return row, files


def record_inspection(row, quiet=False):
    """Accumulate one row per series, newest value wins."""
    os.makedirs(GEO_REC, exist_ok=True)
    path = os.path.join(GEO_REC, "inspected.csv")
    df = pd.DataFrame([row])
    if os.path.exists(path):
        old = pd.read_csv(path)
        df = pd.concat([old[old["accession"] != row["accession"]], df],
                       ignore_index=True)
    df.to_csv(path, index=False)
    if not quiet:
        print(f"  Recorded in curation/GEO/inspected.csv")


def cached_reasons(value):
    """Split a stored reason field back into a list.

    Not `str(value or "")`: these rows come back from a CSV, where an empty
    field reads as NaN, and NaN is truthy in Python. `NaN or ""` returns
    NaN, so str() made the literal string "nan" and every survivor picked
    up a phantom rejection reason.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [x for x in str(value).split(" | ") if x and x != "nan"]


def decide(row, target, include_flagged):
    """Run the two expensive gates on a survivor and settle its verdict.

    Returns (decision, hard, soft, info) where decision is DOWNLOAD,
    AMBIGUOUS or SKIP. DOWNLOAD means the samples declared tissue and
    nothing contradicted it; AMBIGUOUS means something could not be
    confirmed, so it is fetched and marked for review.
    """
    gse_id = row["accession"]
    files, _ = list_suppl_files(gse_id)

    h3, s3, i3 = gate_contents(gse_id, files, row.get("n_samples") or 0)
    time.sleep(NCBI_PAUSE)
    h4, s4, i4 = gate_samples(gse_id, target, include_flagged)

    hard = cached_reasons(row.get("hard")) + h3 + h4
    soft = cached_reasons(row.get("soft")) + s3 + s4

    # Tissue confirmed settles the small-sample suspicion, which was only
    # ever a proxy for "probably cultured": GSE244331 has 7 samples that
    # declare "tumor tissue" and belongs in the analysis.
    if i4.get("declared") == "tissue":
        soft = [x for x in soft
                if "samples (<" not in x and "the abstract mentions" not in x]

    if hard:
        return "SKIP", hard, soft, {**i3, **i4}
    if soft:
        return "AMBIGUOUS", hard, soft, {**i3, **i4}
    return "DOWNLOAD", hard, soft, {**i3, **i4}


# ══════════════════════════════════════════════════════════════════════════
#  Download
# ══════════════════════════════════════════════════════════════════════════

def download_suppl(gse_id, raw_dir, files):
    """Fetch each file by its own URL.

    Not a recursive crawl: ftp.ncbi.nlm.nih.gov serves robots.txt with
    "Disallow: /", so wget -r follows nothing and comes back empty. The
    file list is already known, which also lets a failure be named.
    """
    base = f"{FTP_HTTPS}/{gse_id[:-3]}nnn/{gse_id}/suppl/"
    failed = []
    for f in files:
        out = os.path.join(raw_dir, f["file"])
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"      {f['file']} exists ({human_size(os.path.getsize(out))})")
            continue
        print(f"      {f['file']} ({f['size']})")
        rc = subprocess.run(["wget", "-c", "-q", "--show-progress",
                             "-O", out, base + f["file"]]).returncode
        if rc != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
            failed.append(f["file"])
            if os.path.exists(out) and os.path.getsize(out) == 0:
                os.remove(out)
    return failed


def sample_metadata(gse_id, out_dir):
    """Per-sample annotations via GEOparse.

    For sequencing series the GSM records carry no expression table — the
    counts live in the supplementary files — so this collects annotations
    only, which is what those tables are reliable for.
    """
    try:
        import GEOparse
    except ImportError:
        print("    GEOparse not installed — skipping sample metadata")
        return
    tmp = os.path.join(out_dir, "_tmp")
    os.makedirs(tmp, exist_ok=True)
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=tmp, silent=True)
        rows = []
        for gsm_id, gsm in gse.gsms.items():
            row = {"gsm_id": gsm_id}
            for key, vals in gsm.metadata.items():
                row[key] = ("; ".join(map(str, vals))
                            if isinstance(vals, list) else str(vals))
            rows.append(row)
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(out_dir, "sample_metadata.tsv"), sep="\t", index=False)
        print(f"    Sample metadata: {len(df)} samples x {len(df.columns)} fields")
    except Exception as e:
        print(f"    ERROR reading series: {e}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def write_summary(out_dir, gse_id, row, info, soft, files):
    """The dataset's card. An ambiguous one leads with how to delete it."""
    rel = os.path.relpath(out_dir, os.path.dirname(DATASET_DIR))
    with open(os.path.join(out_dir, "summary.md"), "w") as f:
        f.write(f"# {gse_id} — {row['title']}\n\n")
        if soft:
            f.write("## NEEDS REVIEW — this dataset is ambiguous\n\n")
            f.write("The filters could not confirm that this series belongs in "
                    "the analysis.\nIt was downloaded so you can look at it, not "
                    "because it passed.\n\n")
            f.write("**Why:**\n\n")
            for x in soft:
                f.write(f"- {x}\n")
            f.write(f"\nCheck the samples on the GEO page below. If it is not "
                    f"needed, delete it:\n\n```bash\nrm -rf {rel}\n```\n\n---\n\n")
        f.write(f"- **Review status:** "
                f"{'AMBIGUOUS' if soft else 'tissue confirmed, passed every gate'}\n")
        f.write(f"- **Inferred cancer:** {row['inferred_cancer']}\n")
        f.write(f"- **Samples:** {row['n_samples']}\n")
        f.write(f"- **Platform:** {row['platform']}\n")
        f.write(f"- **Samples declare:** {info.get('declared', '?')} "
                f"— `{info.get('sample_detail', '')}`\n")
        f.write(f"- **Counts are whole numbers:** "
                f"{info.get('raw_is_integer') or 'not verifiable'}\n")
        f.write(f"- **Gene identifiers:** {info.get('gene_ids') or '?'}\n")
        f.write(f"- **Matrix columns:** {info.get('matrix_columns') or '?'}\n\n")
        f.write("## Supplementary files\n\n")
        for fl in files:
            f.write(f"- `{fl['file']}` ({fl['size']}) — {fl['kind']}\n")
        f.write(f"\n## Source\n{row['url']}\n")


def download_gse(gse_id, out_root, row, info, soft, quiet=False):
    out_dir = os.path.join(out_root, gse_id)
    raw_dir = os.path.join(out_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    files, _ = list_suppl_files(gse_id)
    print("    Sample annotations...")
    sample_metadata(gse_id, out_dir)

    # A <GSE>_RAW.tar only repackages the per-sample files, so when a count
    # matrix is published too the bundle is redundant and usually enormous:
    # GSE162739 ships 7.5 GB of bundle beside a 1 MB Raw_counts.tsv.gz.
    wanted = files
    bundles = [f for f in files if f["kind"] == "bundle"]
    if bundles and any(f["kind"] == "raw" for f in files):
        saved = sum(parse_size(f["size"]) for f in bundles)
        print(f"    Skipping {len(bundles)} redundant bundle(s), "
              f"{human_size(saved)}: a count matrix is published separately")
        wanted = [f for f in files if f["kind"] != "bundle"]

    print(f"    Supplementary files ({len(wanted)})...")
    failed = download_suppl(gse_id, raw_dir, wanted)
    got = [f for f in os.listdir(raw_dir) if not f.startswith(".")]
    total = sum(os.path.getsize(os.path.join(raw_dir, f)) for f in got)
    print(f"    On disk: {len(got)} files ({human_size(total)})")
    if failed:
        print(f"    FAILED: {', '.join(failed)} — fetch by hand from {row['url']}")

    if files:
        pd.DataFrame(files).to_csv(os.path.join(out_dir, "suppl_files.csv"),
                                   index=False)
    write_summary(out_dir, gse_id, row, info, soft, files)
    print(f"    summary.md written"
          f"{' (NEEDS REVIEW)' if soft else ''}")
    return out_dir


# ══════════════════════════════════════════════════════════════════════════
#  Pipeline
# ══════════════════════════════════════════════════════════════════════════

def run_pipeline(target, max_download=None, include_flagged=False,
                 min_samples=10):
    cancers = (CANCER_TYPES if target == "ALL" else
               {k: v for k, v in CANCER_TYPES.items() if v["short"] == target})
    if not cancers:
        print(f"Unknown cancer: {target}. Valid: {', '.join(CANCER_SHORTS)}, all")
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"  CLAUD-IA WP1 — GEO pipeline: "
          f"{'all 6 GI cancers' if target == 'ALL' else target}")
    print(f"{'=' * 60}")

    # ── [1/5] Search ────────────────────────────────────────────────────
    print("\n  [1/5] Searching GEO...")
    found = []
    for project_id, info in cancers.items():
        print(f"    {info['label']}")
        for meta in search_candidates(info["keywords"]):
            found.append({"accession": meta["accession"],
                          "search_cancer": info["short"],
                          "n_samples": meta.get("n_samples", ""),
                          "date": meta.get("pdat", ""),
                          "title": str(meta.get("title", ""))[:140],
                          "url": f"{ACC_CGI}?acc={meta['accession']}"})
        time.sleep(1)
    if not found:
        print("    No results.")
        return

    search_df = pd.DataFrame(found).drop_duplicates(subset="accession")
    for short, part in search_df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_rec(short), "search_results.csv"),
                    index=False)
    cancer_of = dict(zip(search_df["accession"], search_df["search_cancer"]))
    print(f"    {len(search_df)} unique candidates")

    # ── [2/5] Gates 1-2 over every candidate ────────────────────────────
    print("\n  [2/5] Gates 1-2: abstract and file names...")
    insp_csv = os.path.join(GEO_REC, "inspected.csv")
    cached = set()
    if os.path.exists(insp_csv):
        cached = set(pd.read_csv(insp_csv)["accession"])
    todo = [a for a in search_df["accession"] if a not in cached]
    print(f"    {len(search_df) - len(todo)} already inspected, {len(todo)} to go")
    for i, acc in enumerate(todo, start=1):
        time.sleep(NCBI_PAUSE)
        try:
            r, _ = inspect_gse(acc, cancer_of.get(acc, target), min_samples,
                               quiet=True)
            verdict = ("rejected" if r["hard"] else
                       "unclear" if r["soft"] else "passes")
            print(f"    [{i}/{len(todo)}] {acc:12} {verdict:9} "
                  f"{(r['hard'] or r['soft'] or '')[:72]}")
        except Exception as e:
            print(f"    [{i}/{len(todo)}] {acc:12} ERROR {e}")

    insp = pd.read_csv(insp_csv)
    insp = insp[insp["accession"].isin(set(search_df["accession"]))]

    # ── [3/5] Gates 3-4 on the survivors ────────────────────────────────
    print("\n  [3/5] Gates 3-4: file contents and declared sample type...")
    survivors = insp[insp["hard"].isna() | (insp["hard"] == "")]
    print(f"    {len(insp)} inspected, {len(survivors)} reached gate 3")

    decisions, payload = [], {}
    for _, row in survivors.iterrows():
        acc = row["accession"]
        time.sleep(NCBI_PAUSE)
        try:
            decision, hard, soft, info = decide(
                row, cancer_of.get(acc, target), include_flagged)
        except Exception as e:
            print(f"      {acc:12} ERROR {e}")
            continue
        payload[acc] = (row, info, soft)
        decisions.append({"accession": acc, "decision": decision,
                          "cancer": row["inferred_cancer"],
                          "n_samples": row["n_samples"],
                          "size": human_size(row["download_bytes"]),
                          "declares": info.get("declared", ""),
                          "reason": "; ".join(hard or soft) or "tissue confirmed",
                          "ambiguity": "; ".join(soft),
                          "title": str(row["title"])[:60]})
        mark = {"DOWNLOAD": "KEEP", "AMBIGUOUS": " ?  ", "SKIP": "skip"}[decision]
        print(f"      [{mark}] {acc:12} {decisions[-1]['reason'][:78]}")

    # Candidates rejected at gates 1-2 go in the record too
    for _, row in insp[~insp["accession"].isin(survivors["accession"])].iterrows():
        decisions.append({"accession": row["accession"], "decision": "SKIP",
                          "cancer": row["inferred_cancer"],
                          "n_samples": row["n_samples"],
                          "size": human_size(row["download_bytes"]),
                          "declares": "", "reason": row["hard"],
                          "ambiguity": "", "title": str(row["title"])[:60]})

    dec_df = pd.DataFrame(decisions)
    dec_df["search_cancer"] = dec_df["accession"].map(cancer_of)
    for short, part in dec_df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_rec(short), "decisions.csv"), index=False)

    print("\n    Why candidates were rejected:")
    skipped = dec_df[dec_df["decision"] == "SKIP"]
    for reason, n in (skipped["reason"].fillna("")
                      .str.split(";").str[0].value_counts().head(15).items()):
        print(f"      {n:4d}  {str(reason)[:78]}")

    selected = dec_df[dec_df["decision"].isin(["DOWNLOAD", "AMBIGUOUS"])]
    amb = dec_df[dec_df["decision"] == "AMBIGUOUS"]
    if max_download and len(selected) > max_download:
        selected = selected.head(max_download)
        print(f"    capped at {max_download} by --max-download")
    if selected.empty:
        print("\n    Nothing reached a download. Reasons per candidate are in "
              "curation/GEO/<CANCER>/decisions.csv")
        return

    # ── [4/5] Download ──────────────────────────────────────────────────
    print(f"\n  [4/5] Downloading {len(selected)} dataset(s) "
          f"({(dec_df['decision'] == 'DOWNLOAD').sum()} confirmed, "
          f"{len(amb)} ambiguous)...")
    for i, acc in enumerate(selected["accession"], start=1):
        row, info, soft = payload[acc]
        tag = "  [NEEDS REVIEW]" if soft else ""
        print(f"\n    --- [{i}/{len(selected)}] {acc}{tag} ---")
        try:
            download_gse(acc, cancer_dir(cancer_of.get(acc, target)),
                         row, info, soft)
        except Exception as e:
            print(f"    ERROR downloading {acc}: {e}")

    # ── [5/5] Summary ───────────────────────────────────────────────────
    print("\n  [5/5] Summary")
    for short in sorted(set(cancer_of.values())):
        n = len([d for d in os.listdir(cancer_dir(short)) if d.startswith("GSE")])
        print(f"    {short}: {n} datasets on disk")
    print(f"    Funnel: {len(search_df)} found -> {len(survivors)} past gates 1-2 "
          f"-> {(dec_df['decision'] == 'DOWNLOAD').sum()} tissue confirmed "
          f"+ {len(amb)} ambiguous")
    if len(amb):
        print("\n    Marked NEEDS REVIEW:")
        for _, r in amb.iterrows():
            print(f"      {r['accession']:12} {r['reason'][:74]}")
    print("\n  Next: python scripts/summary_report.py --md5")


def run_search(target):
    """Gate 1's input only: the candidate list, without inspecting."""
    cancers = (CANCER_TYPES if target == "ALL" else
               {k: v for k, v in CANCER_TYPES.items() if v["short"] == target})
    if not cancers:
        print(f"Unknown cancer: {target}. Valid: {', '.join(CANCER_SHORTS)}, all")
        sys.exit(1)
    rows = []
    for project_id, info in cancers.items():
        print(f"\n  Searching GEO: {info['label']}")
        for meta in search_candidates(info["keywords"]):
            rows.append({"accession": meta["accession"],
                         "search_cancer": info["short"],
                         "n_samples": meta.get("n_samples", ""),
                         "date": meta.get("pdat", ""),
                         "title": str(meta.get("title", ""))[:140]})
        time.sleep(1)
    if not rows:
        print("\n  No results.")
        return
    df = pd.DataFrame(rows).drop_duplicates(subset="accession")
    for short, part in df.groupby("search_cancer"):
        part.to_csv(os.path.join(cancer_rec(short), "search_results.csv"),
                    index=False)
    print(f"\n  {len(df)} unique candidates, saved to "
          f"curation/GEO/<CANCER>/search_results.csv")
    print(f"  Next: python scripts/download_geo.py {target}")


# ══════════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════════

def usage():
    shorts = "|".join(CANCER_SHORTS)
    print("Usage:")
    print(f"  python scripts/download_geo.py <{shorts}|all> [options]")
    print("        search -> gates 1-2 -> gates 3-4 -> download -> summary")
    print()
    print("  Options:")
    print("    --max-download N   stop after N datasets")
    print("    --include-flagged  fetch the cell-line/organoid series too")
    print("    --min-samples N    suspicion threshold (default 10, 0 disables)")
    print()
    print("  Single steps:")
    print(f"    search <{shorts}|all>     candidate list only")
    print("    inspect <GSE_ID>              one series, in detail")
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    mode = sys.argv[1]
    start_logging()

    def opt(name, default=None, cast=int):
        return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default

    if mode.upper() in CANCER_SHORTS or mode.upper() == "ALL":
        run_pipeline(mode.upper(),
                     max_download=opt("--max-download"),
                     include_flagged="--include-flagged" in sys.argv,
                     min_samples=opt("--min-samples", 10))
    elif mode.lower() == "search" and len(sys.argv) > 2:
        run_search(sys.argv[2].upper())
    elif mode.lower() == "inspect" and len(sys.argv) > 2:
        inspect_gse(sys.argv[2].upper())
    else:
        usage()
