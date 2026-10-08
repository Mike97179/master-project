#!/usr/bin/env python3
"""Download the datasets that the decimal rule rejected by mistake.

download_geo.py treats a decimal value in a "raw counts" file as proof that
the file is normalised, and rejects on the spot. That is right for FPKM, TPM
and CPM, where the sequencing depth has already been divided out and cannot be
recovered. It is wrong for RSEM and salmon.

Those tools split a read that fits several isoforms of the same gene in
proportion to its likelihood, so a gene can legitimately receive 0.7 of a
read. The result keeps the scale of a count and its relation to library
depth, which is what DESeq2 needs; rounding it is the standard route, and what
tximport does. The rule could not tell the two apart.

A manual review of the 54 series rejected for that reason alone found 36 with
declared human tissue and recoverable counts: 2,527 samples. This downloads
exactly those, listed in curation/GEO/decimal_rejections.csv.

It does not touch the pipeline, its cache, or any existing decision. Each
series lands in its own cancer folder beside the rest, with a summary.md that
states why it is here and what has to be done to its counts before use.

    python scripts/recover_decimals.py              # all of them
    python scripts/recover_decimals.py PAAD COAD    # only these cancers
    python scripts/recover_decimals.py --dry-run    # just list the plan
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import download_geo as dg
from download_geo import GEO_ROOT, GEO_REC
from run_logger import start_logging

LIST = os.path.join(GEO_REC, "decimal_rejections.csv")

# Files whose columns are counts and FPKM side by side. The decimals the
# detector found were in the FPKM columns; the count columns are whole. These
# need column selection before use, so they are called out by name.
INTERLEAVED = {
    "GSE146889": "keep the columns ending in `_count`, drop the `_rpkm` ones",
    "GSE92945":  "keep the bare sample columns, drop the ones ending in `_fpkm`",
}

# Series whose matrix carries clinical rows or annotation columns ahead of the
# data, so a naive read_csv will not give a count matrix.
STRUCTURED = {
    "GSE254660": "the first rows are clinical annotation prefixed with `#` "
                 "(age_at_diagnosis, alivedeath, recurrence); skip them",
    "GSE208732": "14 annotation columns precede the samples (coordinates, "
                 "HGNC, gene name)",
    "GSE243584": "same layout as GSE208732: annotation columns first",
    "GSE146889": "7 annotation columns precede the samples",
}


def load_plan(cancers):
    if not os.path.exists(LIST):
        sys.exit(f"  Missing {LIST} — see curation/GEO/REVISION_DECIMALES.md")
    d = pd.read_csv(LIST)
    d = d[(d["veredicto"] == "RECUPERABLE") & (d["declara"] == "tissue")]
    if cancers:
        d = d[d["cancer"].isin(cancers)]
    return d.sort_values("n", ascending=False)


def note_for(acc):
    """What has to be done to this dataset's counts before using them."""
    notes = []
    if acc in INTERLEAVED:
        notes.append(f"**Interleaved columns.** {INTERLEAVED[acc]}")
    if acc in STRUCTURED:
        notes.append(f"**Non-standard layout.** {STRUCTURED[acc]}")
    if acc not in INTERLEAVED:
        notes.append("**Round the values before DESeq2.** These are RSEM or "
                     "salmon expected counts: `round()` them, as tximport "
                     "does. Document the rounding in the write-up.")
    return notes


def rewrite_summary(out_dir, acc, n, cancer):
    """Replace write_summary's NEEDS REVIEW boilerplate.

    download_geo writes "the filters could not confirm that this series
    belongs in the analysis... not because it passed". For these 36 that is
    false: they were confirmed by hand. The notes are processing
    instructions, not doubts, and leaving the generic text in place would
    make the record say the opposite of what happened.
    """
    path = os.path.join(out_dir, "summary.md")
    if not os.path.exists(path):
        return
    t = open(path, encoding="utf-8").read()
    if "## NEEDS REVIEW" not in t:
        return
    ini = t.index("## NEEDS REVIEW")
    fin = t.index("---", t.index("rm -rf"))
    notes = "\n".join(f"- {x}" for x in note_for(acc))
    t = t[:ini] + (
        f"## RECOVERED 2026-10-07 — rejected by mistake, confirmed by hand\n\n"
        f"The pipeline rejected this series with *\"the 'raw' counts contain "
        f"decimals, so they are not counts\"*. That rule cannot tell FPKM and "
        f"TPM, where the sequencing depth is gone for good, from the expected "
        f"counts of RSEM and salmon, which split reads between isoforms in "
        f"proportion to their likelihood and keep the scale of a count.\n\n"
        f"The matrix was opened and checked: the values have the magnitude of "
        f"counts, and gate 4 had already confirmed `declares: tissue` for all "
        f"{int(n)} samples. **Treat this as confirmed**, not ambiguous.\n\n"
        f"### Before you use the counts\n\n{notes}\n\n"
        f"Full reasoning and the measurements behind it: "
        f"`curation/GEO/REVISION_DECIMALES.md`\n"
        f"Row for this series: `curation/GEO/decimal_rejections.csv`\n\n"
    ) + t[fin:]
    open(path, "w", encoding="utf-8").write(t)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    fix = "--fix-summaries" in sys.argv
    start_logging()

    plan = load_plan([a.upper() for a in args])
    print("\n" + "=" * 60)
    print("  CLAUD-IA WP1 — recovering the decimal rejections")
    print("=" * 60)
    print(f"\n  {len(plan)} datasets, {int(plan['n'].sum())} samples")
    for c, g in plan.groupby("cancer"):
        print(f"    {c:5} {len(g):2} datasets, {int(g['n'].sum()):4} samples")

    if fix:
        print("\n  Rewriting the summaries of what is already on disk...")
        for _, r in plan.iterrows():
            d = os.path.join(GEO_ROOT, r["cancer"], r["accession"])
            if os.path.isdir(d):
                rewrite_summary(d, r["accession"], r["n"], r["cancer"])
                print(f"    {r['accession']}")
        return

    if dry:
        print("\n  --dry-run: nothing downloaded\n")
        print(plan[["accession", "cancer", "n", "file"]]
              .to_string(index=False, max_colwidth=50))
        return

    insp = pd.read_csv(os.path.join(GEO_REC, "inspected.csv")).set_index("accession")
    done, failed = [], []

    for i, (_, r) in enumerate(plan.iterrows(), start=1):
        acc = r["accession"]
        out_root = os.path.join(GEO_ROOT, r["cancer"])
        print(f"\n  --- [{i}/{len(plan)}] {acc}  ({r['cancer']}, "
              f"{int(r['n'])} samples) ---")
        if os.path.isdir(os.path.join(out_root, acc, "raw")):
            print("    already on disk, skipping")
            done.append(acc)
            continue
        try:
            row = insp.loc[acc].to_dict()
            row["accession"] = acc
            info = {"raw_is_integer": "no (RSEM/salmon expected counts)",
                    "gene_ids": row.get("gene_ids", ""),
                    "matrix_columns": row.get("matrix_columns", "")}
            time.sleep(dg.NCBI_PAUSE)
            dg.download_gse(acc, out_root, row, info, note_for(acc))
            rewrite_summary(os.path.join(out_root, acc), acc,
                            r['n'], r['cancer'])
            done.append(acc)
        except Exception as e:
            print(f"    ERROR {type(e).__name__}: {e}")
            failed.append(acc)

    print("\n" + "=" * 60)
    print(f"  Recovered {len(done)} of {len(plan)}")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
        print("  Rerun to retry: datasets already on disk are skipped")
    print("\n  Next: python scripts/summary_report.py --md5")


if __name__ == "__main__":
    main()
