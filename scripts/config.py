"""
config.py — Shared configuration for CLAUD-IA WP1 data collection.

Centralizes cancer types, database URLs, sample filters, and paths
so every script uses the same definitions.
"""

import os

# ── Project root (one level up from scripts/) ─────────────────────────────
PROJECT_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
DATASET_DIR = os.path.join(PROJECT_ROOT, "dataset")

# ── GI cancer types ───────────────────────────────────────────────────────
# Keys = TCGA project IDs; values = search terms for GEO/PubMed
CANCER_TYPES = {
    "TCGA-COAD": {
        "label": "Colorectal / Colon Cancer",
        "short": "COAD",
        "keywords": ["colon cancer", "colorectal cancer", "rectal cancer"],
    },
    "TCGA-PAAD": {
        "label": "Pancreatic Cancer / PDAC",
        "short": "PAAD",
        "keywords": ["pancreatic cancer", "pancreatic ductal adenocarcinoma", "PDAC"],
    },
    "TCGA-ESCA": {
        "label": "Esophageal Cancer",
        "short": "ESCA",
        "keywords": ["esophageal cancer", "esophageal carcinoma",
                      "gastroesophageal junction cancer", "oesophageal cancer"],
    },
    "TCGA-STAD": {
        "label": "Gastric / Stomach Cancer",
        "short": "STAD",
        "keywords": ["gastric cancer", "stomach cancer"],
    },
    "TCGA-LIHC": {
        "label": "Liver Cancer / HCC",
        "short": "LIHC",
        "keywords": ["liver cancer", "hepatocellular carcinoma", "HCC"],
    },
    "TCGA-CHOL": {
        "label": "Biliary Cancer / Cholangiocarcinoma",
        "short": "CHOL",
        "keywords": ["cholangiocarcinoma", "biliary cancer", "bile duct cancer",
                      "gallbladder cancer"],
    },
}

# ── Cancer inference patterns ─────────────────────────────────────────────
# The keywords above are search terms; these are for recognising which cancer
# a paper is about. They have to be separate: "esophageal adenocarcinoma"
# does not contain the phrase "esophageal cancer", so exact-phrase matching
# labelled a third of the ESCA hits UNKNOWN.
#
# Ordered most specific first, which breaks ties: intrahepatic
# cholangiocarcinoma mentions both bile ducts and liver, and CHOL is the
# more specific answer.
CANCER_ORGAN_PATTERNS = {
    "CHOL": r"cholangio|biliary|bile[-\s]?duct|gall[-\s]?bladder|\bicca\b|\bgbc\b",
    "ESCA": r"esophag|oesophag|barrett|\besca\b|\bescc\b|\beac\b|\bgej\b",
    "PAAD": r"pancrea|\bpdac\b",
    "COAD": r"colorectal|\bcolon\b|\brect(al|um)\b|\bcrc\b",
    "STAD": r"gastric|stomach",
    "LIHC": r"hepatocellular|hepatic|\bliver\b|\bhcc\b",
}

# ── GTEx tissues → WP1 cancer type mapping ────────────────────────────────
GTEX_TISSUE_MAP = {
    "Pancreas":                                "PAAD",
    "Colon - Transverse":                      "COAD",
    "Colon - Sigmoid":                         "COAD",
    "Stomach":                                 "STAD",
    "Esophagus - Mucosa":                      "ESCA",
    "Esophagus - Muscularis":                  "ESCA",
    "Esophagus - Gastroesophageal Junction":   "ESCA",
    "Liver":                                   "LIHC_CHOL",
    "Small Intestine - Terminal Ileum":         "INTESTINE",
}

# ── Sample type filter (Abrar directive) ──────────────────────────────────
# TCGA sample type codes: 01 = Primary Tumor, 11 = Solid Tissue Normal
# We EXCLUDE: cell lines, organoids, xenografts
TCGA_VALID_SAMPLE_CODES = {"01", "11", "06"}
TCGA_SAMPLE_LABELS = {
    "01": "Primary Tumor",
    "11": "Solid Tissue Normal",
    "06": "Metastatic",
}

# Full TCGA sample type vocabulary, so an excluded sample can be named in the
# log instead of being lumped under "cell lines". Codes 02 (recurrent) and 05
# (new primary) ARE human biopsies: they are excluded on purpose, because a
# post-treatment or second primary tumour is not comparable to a primary one,
# not because they fail Abrar's biopsy-only directive.
TCGA_SAMPLE_TYPE_CODES = {
    "01": "Primary Solid Tumor",
    "02": "Recurrent Solid Tumor",
    "03": "Primary Blood Derived Cancer - Peripheral Blood",
    "05": "Additional - New Primary",
    "06": "Metastatic",
    "07": "Additional Metastatic",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
    "12": "Buccal Cell Normal",
    "14": "Bone Marrow Normal",
    "20": "Control Analyte",
    "40": "Recurrent Blood Derived Cancer - Peripheral Blood",
    "50": "Cell Lines",
    "60": "Primary Xenograft Tissue",
    "61": "Cell Line Derived Xenograft Tissue",
}

# ── GDC API endpoints ────────────────────────────────────────────────────
GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"

# ── TCGA RNA-seq filter (bulk RNA-seq only, per Abrar) ───────────────────
TCGA_RNASEQ_FILTER = {
    "data_category": "Transcriptome Profiling",
    "data_type": "Gene Expression Quantification",
    "workflow_type": "STAR - Counts",
    "access": "open",
}
