# Single cgMLST target alleles identify the linezolid-resistant *Staphylococcus capitis* L-clone

Data and analysis code for a study asking whether a single core-genome MLST
(cgMLST) target allele can identify the linezolid-resistant *S. capitis*
lineage (Cluster L), and whether such a marker keeps its performance in genome
collections other than the one it was derived from.

Manuscript in preparation. Citation details will be added on acceptance.

## What the analysis does

1. Every allele at each of the 1,492 targets of the *S. capitis* cgMLST scheme
   is cross-tabulated against Cluster L membership in the 602 labelled genomes
   of a 620-genome public discovery collection.
2. Alleles present in all 19 Cluster L genomes are retained as candidate
   markers (1,346 of 1,492 targets) and ranked by the number of non-lineage
   carriers.
3. Every candidate is then scored in an independent collection of 301 genomes
   deposited from 2023 onward, screened for independence from the discovery
   collection on accession stem, BioSample and whole-profile comparison.
4. *ftsZ* is examined in detail because it belongs to an established Sanger
   sequencing typing scheme and can therefore be read without whole-genome
   sequencing. Its coding sequence is extracted from 1,039 assemblies and the
   resulting haplotypes are compared with the cgMLST allele partition.

## Layout

```
data/        raw inputs (see Provenance below)
scripts/     analysis scripts, in the order given by run_all.sh
results/
  tables/    derived tables, plus the manuscript tables (Tables_1-4_and_S1-S3.xlsx)
  figures/   Figures 1-3 as PNG and PDF, with their legends
run_all.sh   driver that runs the whole analysis in order
```

## Reproducing

```bash
pip install pandas numpy openpyxl matplotlib biopython
bash run_all.sh
```

Scripts expect to be run from the repository root. Several take their inputs as
command-line arguments; `run_all.sh` supplies them.

Run with: Python 3.13.18, pandas 3.0.5, numpy 2.5.3, openpyxl 3.1.5,
matplotlib 3.11.2, biopython 1.88. The discovery screen reproduces 1,346
candidate targets of 1,492, 26 of them with no false positive, and *ftsZ*
allele 4 at rank 50 with 11 false positives (specificity 98.11%, PPV 63.3%).

## Provenance

All genome data are public. Nothing here derives from patient records, and no
clinical isolate is included.

| File | Source |
|---|---|
| `data/Scapitis_620_cgMLST_2026.xlsx` | cgMLST allele calls for the 620 discovery genomes, produced with MBioSEQ Ridom Typer v12.0.5 (Ridom GmbH, Münster, Germany) under the *S. capitis* cgMLST scheme of Wang et al., Nat Commun 2022;13:4254 |
| `data/Scapitis_Merged_Standardized.csv` | cgMLST allele calls for the 301 external genomes, same scheme and software |
| `data/wang2022_cgmlst_clusters_602.csv` | published cgMLST cluster labels, from the supplementary material of Wang et al., Nat Commun 2022;13:4254 |
| `data/ftsZ_CDS_all.fasta` | *ftsZ* coding sequences extracted de novo from 1,039 assemblies in this study (1,188 bp each) |

Genome assemblies and sequencing reads were obtained from NCBI. The seven-locus
MLST scheme used for sequence typing is that of Wang et al.,
BMC Microbiol 2025;25:645.

## Relationship to the companion study

The 620 discovery genomes are also analysed in a separate submitted manuscript
by the same group comparing multilocus sequence typing schemes for
*S. capitis* (repository:
[scapitis-mlst-validation](https://github.com/MoanShane/scapitis-mlst-validation)).
The two studies share that public genome set and nothing else.

Deliberately not in this repository, because it belongs to the companion study:
the six-locus allele catalogue, the provisional sequence types, the clinical
isolate series, the six-gene phylogeny, and the reconciliation of allele
numbering between schemes.

## License

Split by content type:

- **Code** (`scripts/`, `run_all.sh`) — MIT License. See `LICENSE`.
- **Data and derived results** (`data/`, `results/`) — Creative Commons
  Attribution 4.0 International (CC BY 4.0). See `LICENSE-DATA`.

The allele definitions of the cgMLST and seven-locus schemes are curated and
distributed by Ridom within MBioSEQ Ridom Typer and are available to licensed
users of the software; they are not redistributed here. Only the allele calls
for the genomes analysed in this study are included.
