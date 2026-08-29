#!/bin/bash
#SBATCH -J g3_proj
#SBATCH -p med-n16-64g
#SBATCH -a 0-99
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 02:00:00
#SBATCH -o logs/stage3_%a.out
#SBATCH -e logs/stage3_%a.err
set -euo pipefail
export P=/fslustre/qhs/ext_fadra_numrah_mayo_edu/BCM_hackathon/chimp_human_correspondence
export PATH=$P/bin:$PATH
export TMPDIR=$P/tmp
cd $P
python3 scripts/stage3_project.py 03_stage2/jobs.tsv $SLURM_ARRAY_TASK_ID \
  00_ref/human/tx_all.fa 00_ref/chimp/chimpanzee.numeric.fa 04_stage3/chunks
