# **IsoChimp: Population-scale isoform diversity exploration in chimpanzee**
<p align="center">
  <img src="IsoChimp.png" alt="IsoChimp Logo" width="300">
</p>

## **Motivation**

The chimpanzee (*Pan troglodytes*) is a key species for comparative and evolutionary genomics, but its transcriptome annotation is fragmented across providers — NCBI RefSeq, Ensembl, and others — each with its own gene models, identifiers, and inclusion criteria. No systematic comparison of these resources currently exists, there is no way to ask a basic question about any chimpanzee transcript — **how common is it across individuals?** — and there is no transcript-level correspondence to the human annotation on the current chimpanzee assembly.

This project addresses all three gaps, as three goals:

1. **Annotation integration:** Merge the major public chimpanzee annotations into a single unified transcript set, and quantify how much they agree and disagree.  
2. **Population-scale isoform index:** Build a read-level index of transcript evidence from public long-read RNA-seq, so that any transcript can be queried for its population frequency and the samples supporting it.  
3. **Human–chimpanzee correspondence:** Map human transcripts onto the chimpanzee assembly and match them against the unified set, so human and chimpanzee isoforms can be compared directly.

The result is a reference transcript set, a queryable population-frequency index, and a human–chimpanzee transcript correspondence table — together letting other groups assess the prevalence of transcripts they detect in their own chimpanzee samples and relate them to human isoforms.

Goals 1 and 2 are independent and run in parallel while goal 3 depends on the output of goal 1.

Flowchart

![](flowchart.png)


## **Reference assembly and conventions**

All work is anchored to a single assembly to avoid coordinate mismatches:

| Item | Value |
| ----- | ----- |
| Assembly | NHGRI\_mPanTro3-v2.1\_pri |
| RefSeq accession | GCF\_028858775.2 |
| GenBank accession | GCA\_028858775.3 |

* Genome FASTA: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz`  
* Sequence-name alias table (RefSeq `NC_*` / GenBank `CM_*` / UCSC `chr*`): `https://hgdownload.soe.ucsc.edu/hubs/GCF/028/858/775/GCF_028858775.2/GCF_028858775.2.chromAlias.txt`

**Naming convention.** All annotations, genome FASTAs, and alignments are normalised to Ensembl/GENCODE-style sequence names (bare numbers, no `chr` prefix) via the chromAlias table before any comparison or merging step. This includes renaming the chimpanzee genome FASTA itself, since read alignments (goal 2\) and lifted annotations (goal 3\) inherit their sequence names from it. Scaffolds absent from the alias table retain their RefSeq accession; the count of such records is reported so the loss is explicit.

---

## **Goal 1 — Chimpanzee annotation integration**

**Goal:** Produce a unified transcript set from multiple reference annotations, and systematically characterize where those annotations agree and where they diverge.

**Input**

* NCBI RefSeq annotation (release RS\_2026\_05): `https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9598/GCF_028858775.2-RS_2026_05/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.gtf.gz`  
* Ensembl annotation (geneset 2025\_05, on GCA\_028858775.3): `https://ftp.ebi.ac.uk/pub/ensemblorganisms/Pan_troglodytes/GCA_028858775.3/ensembl/geneset/2025_05/genes.gtf.gz`  
* Reference genome FASTA and chromAlias table (above)  
* `isomatch` — https://github.com/zhengxinchang/isomatch

**Steps.**

1. **Download the inputs**  
   1. NCBI RefSeq GTF  
   2. Ensembl GTF  
   3. Chimpanzee reference genome FASTA  
   4. `chromAlias` table  
2. **Standardize chromosome names**  
   1. Normalize chromosome names in both GTFs and the reference genome to the same convention.  
   2. Use **Ensembl-style chromosome names** (e.g. `1`, `2`, `X`, without the `chr` prefix).  
   3. Use the `chromAlias` table when necessary to map RefSeq/GenBank/UCSC chromosome names.  
3. **Merge the two annotations with IsoMatch**  
   1. Run the IsoMatch `merge` function using the two normalized GTFs and the normalized reference genome.  
   2. Keep RefSeq and Ensembl as separate tracked sources.  
   3. The merged GTF should contain provenance information such as `ISOM_SRC` and other IsoMatch attributes.  
4. **Parse the merged GTF and summarize annotation overlap**  
   Calculate at least:  
   1. Total number of transcripts in RefSeq and Ensembl  
   2. Number of **mono-exon** and **multi-exon** transcripts in each annotation  
   3. Number and percentage of transcripts **shared between RefSeq and Ensembl**  
   4. Number of transcripts **unique to RefSeq**  
   5. Number of transcripts **unique to Ensembl**  
   6. If possible, distinguish **exact intron-chain matches** from transcripts with the same splice structure but different transcript boundaries.  
5. **Generate summary outputs**  
   1. A unified chimpanzee annotation GTF with source provenance  
   2. A small summary table of transcript counts and overlap  
   3. An **UpSet plot or Venn diagram** showing RefSeq/Ensembl shared and unique transcripts

**Output**

1. A unified chimpanzee transcript set (GTF) with provenance for every transcript.  
2. Cross-annotation comparison figures — an UpSet plot (or Venn diagram) of transcript overlap between sources, plus a short summary table of per-source counts.

---

## **Goal 2 — Read-level, population-scale isoform index**

**Goal:** Build an index of transcript evidence at the level of individual long reads, aggregated across all publicly available chimpanzee long-read RNA-seq samples, supporting population-frequency queries for arbitrary query transcripts.

**Input**

* Public chimpanzee long-read RNA-seq runs from SRA, as a fixed list of run accessions (BioProject IDs and the accession list are recorded in the project repository; session links to the SRA Run Selector are not stable and are not used)  
* Reference genome, normalised as above  
* `isopedia` — https://github.com/zhengxinchang/isopedia


**Steps.**

1. **Curate chimpanzee long-read RNA-seq datasets**  
   1. Combine SRA run information with sample metadata.  
   2. Identify tissue, individual, platform, and library type.  
   3. Filter out unsuitable runs, including raw PacBio subreads when processed reads are available, datasets without usable transcript reads, and targeted experiments.  
   4. Generate the final SRR list and a QC/audit table.  
2. **Current result: 17 suitable samples identified.**  
3. **Prepare the reference**  
   1. Normalize chromosome names in the chimpanzee genome and GTF to the same convention.  
   2. Build the required `samtools` and `minimap2` indexes.  
4. **Download and align the selected runs**  
   1. Download FASTQ files for the approved SRR accessions.  
   2. Verify downloads with MD5 checksums.  
   3. Align reads to the normalized chimpanzee genome with splice-aware `minimap2`.  
   4. Sort/index BAM files and merge runs belonging to the same BioSample.  
   5. Calculate basic alignment QC, including mapped reads and mapping rate.  
5. **Build the Isopedia index**  
   1. Build a population-scale `isopedia` index from the cleaned BAM files while retaining sample provenance.  
6. **Query transcript set**   
   1. Use the unified chimpanzee GTF from Goal 1 as the query or use ENSEMBL gene.gtfs as an anternative.  
   2. Summarize transcript support, including:  
      1. Supporting reads  
      2. Supporting samples  
      3. Population frequency  
      4. Sample IDs  
7. **Generate summary outputs**  
   1. Curated sample metadata and inclusion/exclusion table  
   2. QCed BAM files  
   3. Chimpanzee `isopedia` population index  
   4. Goal 1 transcripts annotated with population support  
   5. Basic population-frequency and sample-support statistics

**Output**

1. A chimpanzee read-level `isopedia` index supporting transcript queries that return population frequency and the list of supporting samples.  
2. The curated SRA run/sample metadata table used to build it.

---

## Goal 3 — Human–chimpanzee transcript correspondence (gene-anchored)

**Goal**

Establish a **transcript-level correspondence between human and chimpanzee isoforms** on `NHGRI_mPanTro3-v2.1_pri`.

The analysis is **gene-anchored**: first identify the corresponding chimpanzee gene/locus for each human gene, then compare transcript structures within that locus. Splice-structure similarity is the primary criterion, while sequence identity and coverage are used as supporting information.

**Input**

* Human GENCODE comprehensive GTF and GRCh38 reference genome  
* NCBI human–chimpanzee gene orthology assignments  
* Normalized chimpanzee reference genome  
* Unified chimpanzee transcript annotation from **Goal 1**  
* `gffread`  
* `minimap2`  
* `Liftoff` — fallback for genes without an NCBI ortholog  
* `isomatch`  
* `BBSketch` or similar tool — optional conservation screen

**Steps**

1. **Optional sequence-conservation screen**  
   * Extract human transcript cDNAs with `gffread`.  
   * Use `BBSketch` or a similar method to estimate human–chimpanzee sequence conservation.  
   * Use this as a QC/sanity check rather than the formal transcript-matching criterion.  
2. **Build the human–chimpanzee gene projection table**  
   * Map GENCODE genes to Entrez GeneIDs.  
   * Retrieve chimpanzee orthologs from NCBI.  
   * Classify gene relationships as:  
     * 1:1 ortholog  
     * 1 / complex  
     * No assigned ortholog  
   * For genes without an NCBI ortholog, use `Liftoff` as a lower-confidence fallback.  
3. **Project human transcripts onto chimpanzee loci**  
   * Extract human transcript sequences and their corresponding chimpanzee ortholog loci.  
   * Splice-align human transcripts to the assigned chimpanzee locus using `minimap2`.  
   * Convert the alignments into projected chimpanzee-coordinate GTF transcript models.  
   * Record alignment identity, coverage, MAPQ, and other basic QC metrics.  
4. **Compare projected transcripts with Goal 1 annotations**  
   * Compare projected human transcript models against the unified chimpanzee GTF using `isomatch`.  
   * Use **intron-chain/splice-structure similarity as the primary comparison**.  
   * Classify transcripts into categories such as:  
     * Exact intron-chain match  
     * Same splice structure with different transcript ends  
     * Compatible / contained  
     * Structurally divergent  
     * Projected but no annotated chimpanzee isoform  
     * Unalignable / no gene anchor  
5. **Summarize and validate**  
   * Calculate the proportion of human transcripts with conserved chimpanzee counterparts.  
   * Summarize structural-match categories and alignment identity/coverage.  
   * Report transcripts without a counterpart and the reason for failure.

**Output**

1. A **human ↔ chimpanzee transcript correspondence table** containing gene anchor, chimpanzee transcript match, structural category, alignment identity, and coverage.  
2. Human transcript projections in **chimpanzee coordinates (GTF)**.  
3. Summary statistics/figures showing:  
   * Conserved versus divergent transcript structures  
   * Distribution of sequence identity/coverage  
   * Human transcripts without a chimpanzee counterpart, stratified by reason.

