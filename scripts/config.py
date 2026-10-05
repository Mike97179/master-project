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
