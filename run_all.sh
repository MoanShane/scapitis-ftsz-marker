#!/bin/bash
# Reproduce the analysis. Run from the repository root:  bash run_all.sh
set -euo pipefail

X=data/Scapitis_620_cgMLST_2026.xlsx       # discovery cgMLST export (620 genomes)
M=data/Scapitis_Merged_Standardized.csv    # external cgMLST export (301 genomes)
C=data/wang2022_cgmlst_clusters_602.csv    # published cluster labels (602 genomes)
F=data/ftsZ_CDS_all.fasta                  # extracted ftsZ coding sequences
A=results/tables/ftsz_allele_assignments.csv

mkdir -p results/tables results/figures

echo "== 1. discovery screen of all 1,492 targets"
python3 scripts/scan1492.py        "$X" "$C"

echo "== 2. external evaluation of every candidate"
python3 scripts/extscan.py         "$X" "$M" "$C"

echo "== 3. candidate-set summary"
python3 scripts/summary.py

echo "== 4. cgMLST distances, both collections"
python3 scripts/dist_all.py        "$X" "$M"
python3 scripts/extdist.py         "$X" "$M"

echo "== 5. sensitivity, specificity, PPV (Wilson) and Fisher's exact test"
python3 scripts/wilson.py

echo "== 6. external genomes flagged by the 26 discovery-perfect candidates"
python3 scripts/whofp.py           "$X" "$M" "$C"

echo "== 7. lineage composition and within-lineage distances"
python3 scripts/lclone_dist.py
python3 scripts/check_nn.py
python3 scripts/impact_602.py

echo "== 8. ftsZ cgMLST alleles"
python3 scripts/ftsz_alleles.py    "$X" "$M"

echo "== 9. ftsZ coding-sequence haplotypes vs cgMLST alleles"
python3 scripts/final_cds.py       "$F"
python3 scripts/cds_xtab.py        "$F" "$A"
python3 scripts/cds_marker.py
python3 scripts/verify2.py         "$F"
python3 scripts/frag_check.py

echo "== 10. minimum discriminating SNP set"
python3 scripts/minset.py

echo "== 11. figures"
python3 scripts/fig1.py
python3 scripts/fig2.py
python3 scripts/fig3b.py

echo "== 12. tables, legends and submission package"
python3 scripts/legends.py
python3 scripts/table_s1.py        "$X" "$M" "$C"
python3 scripts/package.py

echo
echo "Done. Outputs are in results/."
