#!/bin/bash
#SBATCH -J g3_stage0
#SBATCH -c 8
#SBATCH -p med-n16-64g
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /fslustre/qhs/ext_fadra_numrah_mayo_edu/BCM_hackathon/chimp_human_correspondence/logs/stage0_%j.out
#SBATCH -e /fslustre/qhs/ext_fadra_numrah_mayo_edu/BCM_hackathon/chimp_human_correspondence/logs/stage0_%j.err
set -euo pipefail
export P=/fslustre/qhs/ext_fadra_numrah_mayo_edu/BCM_hackathon/chimp_human_correspondence
export PATH=$P/bin:$PATH
export TMPDIR=$P/tmp
cd $P
CHIMP_GTF=00_ref/chimp/chimp_refseq.numeric.gtf.gz
zcat $CHIMP_GTF | awk -F"\t" '$3=="exon" && $9 ~ /transcript_id "[^"]+"/' > $P/tmp/chimp_for_gffread.gtf
/usr/local/biotools/gffread/0.12.7/gffread $P/tmp/chimp_for_gffread.gtf -g 00_ref/chimp/chimpanzee.numeric.fa -w 01_stage0/chimp_tx.fa
/usr/local/biotools/bbmap/38.90/sketch.sh in=01_stage0/chimp_tx.fa out=01_stage0/chimp_tx.sketch mode=sequence k=31,24 threads=8 -Xmx16g
/usr/local/biotools/bbmap/38.90/comparesketch.sh in=00_ref/human/tx_all.fa ref=01_stage0/chimp_tx.sketch \
  mode=sequence k=31,24 format=3 records=1 minhits=3 threads=8 \
  out=01_stage0/human_vs_chimp_tx.tsv -Xmx16g
/usr/local/biotools/python/3.10.11/bin/python scripts/stage0_summary.py 01_stage0/human_vs_chimp_tx.tsv 00_ref/human/tx_meta.tsv 01_stage0/stage0
