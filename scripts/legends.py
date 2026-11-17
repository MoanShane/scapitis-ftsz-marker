import os, shutil
from docx import Document
from docx.shared import Pt

PK="results/package"
LEG=[("Figure 1.",
  "Selection of the external genome collection. Of 952 Staphylococcus capitis assembly records in the "
  "NCBI Assembly database, 614 were released from 2023 onward; paired GenBank and RefSeq records were "
  "collapsed to 307 genomes, 5 failed the independence or quality filters, and the S. capitis subsp. "
  "capitis type strain DSM 20326 (GCF_040739495.1) was removed on review as the outgroup of the "
  "discovery collection, giving a final set of 301 genomes."),
 ("Figure 2.",
  "cgMLST allelic distance to the Cluster L reference strain LNZR-1 in the two collections. Each point "
  "is one genome, jittered vertically; the horizontal axis is broken between 60 and 385 allelic "
  "differences, and all 921 genomes are plotted. In the discovery collection (n = 620) the lineage "
  "reaches 25 allelic differences and the nearest non-lineage genome sits at 398, a gap of 373. In the "
  "external collection (n = 301) seven-locus ST6 genomes extend to 29 — three of them beyond the "
  "25-allele cutoff derived from the discovery set — and a single ST46 genome sits at 50, before a gap "
  "of 464 to the next genome at 514. Marker shape as well as colour distinguishes the groups, so the "
  "figure is readable in greyscale and under colour-vision deficiency."),
 ("Figure 3.",
  "Discovery against external performance for all 1,346 candidate targets. Each point is one candidate, "
  "placed by the number of non-lineage genomes it flagged in the discovery collection and the number it "
  "flagged in the external collection, with the 26 discovery-perfect candidates and ftsZ marked. The "
  "figure carries two things at once: how closely discovery performance predicts external performance "
  "(Spearman rho = 0.94), and the fact that no candidate combined complete external sensitivity with "
  "zero external false positives. Both axes use a symmetric log scale so that zero can be shown.")]

d=Document()
st=d.styles["Normal"]; st.font.name="Calibri"; st.font.size=Pt(11)
h=d.add_paragraph(); r=h.add_run("Figure legends"); r.bold=True; r.font.size=Pt(13)
sub=d.add_paragraph()
r2=sub.add_run("Single cgMLST Target Alleles Identify the Linezolid-Resistant Staphylococcus capitis "
               "L-Clone: Population-Dependent Predictive Value  —  manuscript v2.0")
r2.italic=True; r2.font.size=Pt(9.5)
for lab,txt in LEG:
    p=d.add_paragraph(); p.paragraph_format.space_after=Pt(10)
    rb=p.add_run(lab+" "); rb.bold=True
    p.add_run(txt)
d.save(f"{PK}/Figure_legends.docx")

with open(f"{PK}/Figure_legends.md","w") as fh:
    fh.write("# Figure legends\n\n*Manuscript v2.0*\n\n")
    for lab,txt in LEG: fh.write(f"**{lab}** {txt}\n\n")

for n in ("Figure1_selection_flow","Figure2_distance_boundary","Figure3_candidate_transfer"):
    for e in ("png","pdf"):
        s=f"results/tables/{n}.{e}"
        if os.path.exists(s): shutil.copy(s, f"{PK}/figures/{n}.{e}")

with open(f"{PK}/README.txt","w") as fh:
    fh.write(
"JMII submission package — manuscript v2.0\n"
"Single cgMLST Target Alleles Identify the Linezolid-Resistant Staphylococcus capitis\n"
"L-Clone: Population-Dependent Predictive Value\n"
"Generated 2026-11-10\n\n"
"Tables_1-4_and_S1-S3.xlsx   One sheet per table.\n"
"   Table 1   Specificity of the 1,346 candidate targets in the discovery collection\n"
"   Table 2   External false positives by discovery specificity band\n"
"   Table 3   ftsZ marker performance, both collections and both gold standards\n"
"   Table 4   Marker-allele carriage among seven-locus ST3 genomes\n"
"   Table S1  Selection of the external genome collection\n"
"   Table S2  External genomes profile-identical to a discovery genome\n"
"   Table S3  The 26 discovery-perfect targets and their external performance\n\n"
"Figure_legends.docx / .md   Legends for Figures 1-3.\n\n"
"figures/   PNG at 400 dpi for review, PDF (vector) for production.\n"
"   Figure1_selection_flow        genome selection, 952 -> 301\n"
"   Figure2_distance_boundary     cgMLST distance to LNZR-1, both collections\n"
"   Figure3_candidate_transfer    discovery vs external performance of all candidates\n\n"
"Not included: the manuscript text itself, which lives in the Claude doc.\n"
"Still outstanding: author list and contributions, funding, competing interests.\n"
"No IRB approval is required; the study uses only public genome data.\n")
print("package contents:")
for root,_,files in os.walk(PK):
    for f in sorted(files): print("  ", os.path.relpath(os.path.join(root,f), PK))
