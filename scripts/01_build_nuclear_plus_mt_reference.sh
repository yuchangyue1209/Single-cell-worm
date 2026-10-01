#!/usr/bin/env bash
# Required: NUCLEAR_FASTA, NUCLEAR_GTF, MT_FASTA, MT_GTF
set -euo pipefail
CELLRANGER_DIR=${CELLRANGER_DIR:-/work/cyu/scRNA/cellranger-10.0.0}
REFERENCE_ROOT=${REFERENCE_ROOT:-/work/cyu/scRNA/reference_test}
REFERENCE_NAME=${REFERENCE_NAME:-ssol_final_nuclear_plus_mt}
BUILD_DIR=${BUILD_DIR:-/work/cyu/scRNA/reference_test/source_nuclear_plus_mt}
: "${NUCLEAR_FASTA:?Set NUCLEAR_FASTA}"
: "${NUCLEAR_GTF:?Set NUCLEAR_GTF}"
: "${MT_FASTA:?Set MT_FASTA}"
: "${MT_GTF:?Set MT_GTF}"
for f in "${NUCLEAR_FASTA}" "${NUCLEAR_GTF}" "${MT_FASTA}" "${MT_GTF}"; do [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 1; }; done
mkdir -p "${BUILD_DIR}" "${REFERENCE_ROOT}"
FA="${BUILD_DIR}/Ssolidus_nuclear_plus_mt.fa"; GTF="${BUILD_DIR}/Ssolidus_nuclear_plus_mt.gtf"
awk 'FNR==1 && NR!=1 {print ""} {print}' "${NUCLEAR_FASTA}" "${MT_FASTA}" > "${FA}.tmp" && mv "${FA}.tmp" "${FA}"
awk 'NF==0 || /^#/ || NF>=9' "${NUCLEAR_GTF}" "${MT_GTF}" > "${GTF}.tmp" && mv "${GTF}.tmp" "${GTF}"
export PATH="${CELLRANGER_DIR}:${PATH}"
[[ ! -e "${REFERENCE_ROOT}/${REFERENCE_NAME}" ]] || { echo "Reference exists" >&2; exit 1; }
cd "${REFERENCE_ROOT}"
cellranger mkref --genome="${REFERENCE_NAME}" --fasta="${FA}" --genes="${GTF}"
sha256sum "${FA}" "${GTF}" > "${BUILD_DIR}/SHA256SUMS"
