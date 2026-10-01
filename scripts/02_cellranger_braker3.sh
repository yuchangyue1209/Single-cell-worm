#!/usr/bin/env bash
# Usage: bash 02_cellranger_braker3.sh Big|Breeding
set -euo pipefail
SAMPLE=${1:?Usage: $0 Big|Breeding}
case "${SAMPLE}" in Big|Breeding) ;; *) echo "Sample must be Big or Breeding" >&2; exit 1 ;; esac
CELLRANGER_DIR=${CELLRANGER_DIR:-/work/cyu/scRNA/cellranger-10.0.0}
REFERENCE=${REFERENCE:-/work/cyu/scRNA/reference_test/ssol_final_nuclear_plus_mt}
FASTQ_DIR=${FASTQ_DIR:-/mnt/spareHD_2/scRNA}
OUTPUT_ROOT=${OUTPUT_ROOT:-/work/cyu/scRNA}
CORES=${CORES:-16}; MEM_GB=${MEM_GB:-90}
RUN_ID=${RUN_ID:-${SAMPLE}_braker_tsa_utr_tagseq_plus_187_mt_2026}
export PATH="${CELLRANGER_DIR}:${PATH}"
[[ -d "${REFERENCE}" ]] || { echo "Missing reference" >&2; exit 1; }
shopt -s nullglob; r1=("${FASTQ_DIR}/${SAMPLE}"*_R1_001.fastq.gz); r2=("${FASTQ_DIR}/${SAMPLE}"*_R2_001.fastq.gz); shopt -u nullglob
[[ ${#r1[@]} -gt 0 && ${#r1[@]} -eq ${#r2[@]} ]] || { echo "Missing/unmatched FASTQs" >&2; exit 1; }
[[ ! -e "${OUTPUT_ROOT}/${RUN_ID}" ]] || { echo "Output exists" >&2; exit 1; }
cd "${OUTPUT_ROOT}"
cellranger count --id="${RUN_ID}" --transcriptome="${REFERENCE}" --fastqs="${FASTQ_DIR}" --sample="${SAMPLE}" --chemistry=auto --include-introns=true --create-bam=false --localcores="${CORES}" --localmem="${MEM_GB}"
