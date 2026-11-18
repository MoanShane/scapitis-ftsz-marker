import os, shutil, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

PK="results/package"; shutil.rmtree(PK, ignore_errors=True)
os.makedirs(f"{PK}/figures", exist_ok=True)

T1=pd.DataFrame({"Non-lineage carriers (of 583)":["0","1–11","12–50","51–200","more than 200"],
                 "Candidate targets":[26,58,511,64,687],
                 "Share of the 1,346":["1.9%","4.3%","38.0%","4.8%","51.0%"]})
T2=pd.DataFrame({"Non-lineage carriers in discovery":["0","1–11","12–50","51–200","more than 200"],
                 "Candidates":[26,58,511,64,687],
                 "With no more than 2 external false positives":["26 (100%)","40 (69.0%)","2 (0.4%)","0","0"]})
T3=pd.DataFrame({"Measure":["Lineage genomes","Non-lineage genomes","Marker carriers","False positives",
                            "Sensitivity","Specificity","Positive predictive value"],
  "Discovery, n = 602":[19,583,30,11,"100.0 (83.2–100.0)","98.1 (96.7–98.9)","63.3 (45.5–78.1)"],
  "ST6, n = 301":[22,279,24,2,"100.0 (85.1–100.0)","99.3 (97.4–99.8)","91.7 (74.2–97.7)"],
  "≤ 50 alleles, n = 301":[23,278,24,1,"100.0 (85.7–100.0)","99.6 (98.0–99.9)","95.8 (79.8–99.3)"],
  "ST6, n = 295":[22,273,24,2,"100.0 (85.1–100.0)","99.3 (97.4–99.8)","91.7 (74.2–97.7)"],
  "≤ 50 alleles, n = 295":[23,272,24,1,"100.0 (85.7–100.0)","99.6 (97.9–99.9)","95.8 (79.8–99.3)"]})
T4=pd.DataFrame({"Collection":["Discovery","External"],"ST3 genomes":[38,27],
                 "Carrying the marker allele":["12 (31.6%)","1 (3.7%)"],
                 "Carrying the common ST3 allele":[21,22]})
S1=pd.read_csv("results/tables/table_s1_accessions.csv", dtype=str).fillna("")
FIG1=pd.DataFrame({"Stage":["NCBI census, all releases","Released 2023 onward",
                          "One record per genome, RefSeq copy kept",
                          "Independence and quality filters passed",
                          "Discovery outgroup removed on review"],
                 "Genomes":[952,614,307,302,301],
                 "Removed at this stage":["—","338 earlier releases","307 paired GenBank duplicates","5",
                                          "1 (GCF_040739495.1, S. capitis subsp. capitis type strain DSM 20326)"]})
S2=pd.DataFrame({"External genome":["GCF_040736935.1","GCF_050044135.1","GCF_054792995.1",
                                    "GCF_040736585.2","GCF_050042105.1","GCF_050042725.1"],
                 "Identical to":["SRR7939128","SRR8786603","DRR015840","SRR7939170",
                                 "GCA_015645265.1","GCA_015645205.1"],
                 "Why the accession filters missed it":["SRA run, no assembly accession to compare"]*4+
                                 ["Assembly, but a different accession stem"]*2})
f=pd.read_csv("results/tables/cgmlst_target_full.csv")
S3=f[f.disc_FP==0][["target","allele","disc_FP","ext_TP","ext_FN","ext_FP","ext_sens"]].copy()
S3.columns=["cgMLST target","Allele","Discovery false positives","External true positives",
            "External false negatives","External false positives (ST6 standard)","External sensitivity (%)"]
S3=S3.sort_values(["External false negatives","cgMLST target"]).reset_index(drop=True)
S3.insert(0,"#",range(1,len(S3)+1))

TITLES={"Table 1":"Table 1. Specificity of the 1,346 candidate targets in the discovery collection",
 "Table 2":"Table 2. External false positives by discovery specificity band",
 "Table 3":"Table 3. ftsZ marker performance across both collections and both gold standards. Values are percentages with Wilson score 95% confidence intervals.",
 "Table 4":"Table 4. Marker-allele carriage among seven-locus ST3 genomes in the two collections",
 "Table S1":"Supplementary Table S1. GenBank and RefSeq accession numbers for the 620 discovery genomes and the 301 external genomes, with the published cgMLST cluster label where one exists and the seven-locus MLST sequence type.",
 "Figure 1 source data":"Figure 1 source data. Selection of the external genome collection. Quality thresholds: total length 2.2–2.9 Mb; contig N50 ≥ 20 kb; ≤ 300 contigs; CheckM completeness ≥ 90%; CheckM contamination ≤ 5%; NCBI ANI ≥ 95% to S. capitis.",
 "Table S2":"Supplementary Table S2. External genomes found profile-identical to a discovery genome across the 1,492 shared cgMLST targets",
 "Table S3":"Supplementary Table S3. The 26 candidate targets with no false positive in the discovery collection, and their performance in the external collection"}
SHEETS=[("Table 1",T1),("Table 2",T2),("Table 3",T3),("Table 4",T4),
        ("Table S1",S1),("Table S2",S2),("Table S3",S3),
        ("Figure 1 source data",FIG1)]

wb=Workbook(); wb.remove(wb.active)
hdr_f=Font(bold=True, color="FFFFFF"); hdr_fill=PatternFill("solid", fgColor="4F6D8C")
ttl_f=Font(bold=True, size=11)
for name,df in SHEETS:
    ws=wb.create_sheet(name)
    ws.cell(1,1,TITLES[name]).font=ttl_f
    ws.cell(1,1).alignment=Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=1,start_column=1,end_row=1,end_column=max(2,len(df.columns)))
    ws.row_dimensions[1].height=34
    for j,c in enumerate(df.columns,1):
        cell=ws.cell(3,j,str(c)); cell.font=hdr_f; cell.fill=hdr_fill
        cell.alignment=Alignment(wrap_text=True, vertical="center", horizontal="center")
    ws.row_dimensions[3].height=30
    for i,row in enumerate(df.itertuples(index=False),4):
        for j,v in enumerate(row,1):
            ws.cell(i,j,v).alignment=Alignment(vertical="center",
                horizontal="left" if j==1 else "center", wrap_text=(j==len(df.columns)))
    for j,c in enumerate(df.columns,1):
        w=max(len(str(c)), *(len(str(v)) for v in df.iloc[:,j-1])) + 3
        ws.column_dimensions[get_column_letter(j)].width=min(max(w,12), 52)
    ws.freeze_panes="A4"
wb.save(f"{PK}/Tables_1-4_and_S1-S3.xlsx")
print("xlsx written:", [n for n,_ in SHEETS])
print("S3 rows:", len(S3))
