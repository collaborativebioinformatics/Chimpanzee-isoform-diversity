# Goal 3 methods summary

Tools: minimap2 2.28 (splice:hq), isomatch (github.com/zhengxinchang/isomatch), gffread 0.12.7,
samtools 1.20, bedtools 2.27.1, BBTools 38.90 (sketch.sh/comparesketch.sh), python3.10
(pandas, pysam, matplotlib, seaborn).

Data: human GENCODE v50 (annotation GTF, GRCh38 primary assembly, EntrezGene metadata);
chimpanzee NHGRI_mPanTro3-v2.1_pri (GCF_028858775.2) genome and RefSeq annotation RS_2026_05;
NCBI gene_orthologs and Homo_sapiens.gene_info (downloaded 2026-08-28).

Pipeline (scripts in ../scripts, full walkthrough in goal3_full_execution_guide.md):
1. Stage 0 sequence screen: gffread cDNA extraction; BBSketch per-sequence sketches (k=31,24) of
   135,577 chimp RefSeq cDNAs; comparesketch best hit per human transcript (stage0_summary.py).
2. Stage 1 gene anchoring: GENCODE->Entrez via GENCODE metadata + NCBI Ensembl xrefs;
   human->chimp GeneID via NCBI gene_orthologs (17,903 pairs, strictly 1:1); coordinates from the
   RefSeq GTF (chromosomes renamed to assembly names by sequence length; stage1_projection.py).
3. Stage 2/3 projection: per anchored gene, human transcript cDNAs splice-aligned to the ortholog
   locus +/-10 kb (minimap2 -ax splice:hq -uf, locus-restricted); alignments converted to exon
   models in chimp coordinates; tiers A/B/C/D by coverage/identity/MAPQ (stage3_project.py, Slurm array).
4. Stage 4 classification: isomatch classify of projected models vs chimp RefSeq; SQANTI-style
   categories mapped to match classes; novel_in/not_in_catalog reported as "not annotated
   (known/novel junctions)", unresolved pending read evidence (stage4_merge.py).
5. Stage 5 summaries and figures (stage5_plots.py).

Deviations from the written plan: comparison reference is chimp RefSeq (Goal 1 unified set not yet
available; swap ANNOT_GTF and rerun Stages 1.2-1.3, 4, 5); Liftoff fallback and Goal 2 read-support
query not run; CDS-only identity not computed. NCBI orthology yielded no 1:many pairs, so the
ambiguous class is empty.
