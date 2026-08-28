#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   ./align_isoseq.sh input.fastq.gz output.sam

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fa/fastq[.gz]> <output.sam>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

GENOME="/path/to/reference/genome.fa"
THREADS=16

minimap2 \
    -t "${THREADS}" \
    -ax splice:hq \
    -uf \
    --secondary=no \
    "${GENOME}" \
    "${INPUT}" \
    > "${OUTPUT}"
