"""Supplementary Table S1: accession list for both collections.

Usage:
    python3 scripts/table_s1.py <620_cgMLST.xlsx> <external_merged.csv> <clusters_602.csv>

Writes results/tables/table_s1_accessions.csv, which package.py turns into the
"Table S1" sheet of Tables_1-3_and_S1-S3.xlsx.

The external collection is selected exactly as extscan.py selects it, so the
table lists the same 301 genomes the analysis used.
"""
import csv
import re
import sys

import openpyxl

XLSX, EXT, CLU = sys.argv[1], sys.argv[2], sys.argv[3]

BATCH = "Scapitis_ExtVal_304thirdBatch"   # as in extscan.py
OUTGROUP_ACC = "GCF_040739495.1"          # S. capitis subsp. capitis type strain DSM 20326
OUTGROUP_SID = "DSM20326"                 # how the discovery export names it
ANCHOR = "GCA_000712995.1"                # LNZR-1, Cluster L reference anchor
NO_LABEL = ("Not assigned", "Outgroup", "nan")   # as in extscan.py


def canon(v):
    if v is None:
        return ""
    if isinstance(v, float) and v != v:              # NaN
        return ""
    s = str(v).strip()
    if s.lower() in ("", "nan", "none", "-", "?"):
        return ""
    m = re.fullmatch(r"(\d+)\.0+", s)                # 3.0 -> 3
    return m.group(1) if m else s


def norm(h):
    """Collapse whitespace and underscores so both exports spell headers alike."""
    return " ".join(str(h).replace("_", " ").split())


def exact(header, wanted):
    for i, h in enumerate(header):
        if norm(h) == wanted:
            return i
    raise SystemExit("no column named %r; saw %s" % (wanted, [norm(h) for h in header[:8]]))


# ---- published cluster labels, country and year (Wang et al. [8]) ----------
labels = {}
with open(CLU, newline="") as fh:
    for row in csv.DictReader(fh):
        labels[canon(row["Accession"])] = (
            canon(row["Wang_Cluster"]),
            canon(row["Country"]),
            canon(row["Year"]),
        )

rows = []

# ---- discovery collection (620) -------------------------------------------
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
ws = wb[wb.sheetnames[0]]
it = ws.iter_rows(values_only=True)
hdr = list(next(it))
c_id = exact(hdr, "Sample ID (Sample)")
c_st = exact(hdr, "ST (S. capitis MLST)")

for r in it:
    sid = canon(r[c_id])
    if not sid:
        continue
    cluster, country, year = labels.get(sid, ("", "", ""))
    if cluster in NO_LABEL:
        cluster = ""
    acc, note = sid, ""
    if sid == OUTGROUP_SID:
        acc = OUTGROUP_ACC
        note = "Outgroup, type strain DSM 20326"
    elif sid == ANCHOR:
        note = "Cluster L reference anchor, strain LNZR-1"
    if not cluster and not note:
        note = "No published cluster label"
    elif not cluster and sid == ANCHOR:
        note = "Cluster L reference anchor, strain LNZR-1; no published cluster label"
    rows.append(["Discovery", acc, cluster, canon(r[c_st]), country, year, note])

# ---- external collection (301) --------------------------------------------
with open(EXT, newline="") as fh:
    rd = csv.reader(fh)
    hdr = next(rd)
    e_id = exact(hdr, "Sample ID (Sample)")
    e_st = exact(hdr, "ST (S. capitis MLST)")
    e_bt = exact(hdr, "SourceBatch")
    for r in rd:
        sid = canon(r[e_id])
        if not sid or canon(r[e_bt]) != BATCH or sid == OUTGROUP_ACC:
            continue
        rows.append(["External", sid, "", canon(r[e_st]), "", "", ""])

rows.sort(key=lambda x: (x[0] != "Discovery", x[1]))

COLS = ["Collection", "Accession", "Published cgMLST cluster [8]",
        "Seven-locus MLST ST", "Country", "Year", "Note"]
with open("results/tables/table_s1_accessions.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(COLS)
    w.writerows(rows)

disc = [r for r in rows if r[0] == "Discovery"]
ext = [r for r in rows if r[0] == "External"]
assert len(set(r[1] for r in ext)) == len(ext), "duplicate external accession"
print("Table S1 written: %d rows (%d discovery, %d external)" % (len(rows), len(disc), len(ext)))
print("  discovery with a published cluster label:", sum(1 for r in disc if r[2]))
print("  discovery without one:", sum(1 for r in disc if not r[2]))
print("  discovery Cluster L:", sum(1 for r in disc if r[2] == "L"))
print("  external ST6:", sum(1 for r in ext if r[3] == "6"))
odd = [r[1] for r in disc if not re.match(r"^(DRR|ERR|SRR|GCA_|GCF_)", r[1])]
print("  discovery IDs that are not accession-shaped:", odd or "none")
