#!/usr/bin/env bash
set -euo pipefail

# Reproduce the final S. solidus single-cell analysis after Cell Ranger.
# Existing result directories are retained; Seurat output files are overwritten
# only when the corresponding R script is intentionally rerun.

PROJECT=/work/cyu/scRNA
SCRIPTS=${PROJECT}/downstream/scripts
RESULTS=${PROJECT}/downstream/results

BIG_MATRIX=${PROJECT}/Big_braker_tsa_utr_tagseq_plus_187_mt_2026/outs/filtered_feature_bc_matrix
BREEDING_MATRIX=${PROJECT}/Breeding_braker_tsa_utr_tagseq_plus_187_mt_2026/outs/filtered_feature_bc_matrix
FUNCTION_MASTER=/work/cyu/annotation/ssol-annotation/results/final_filtered_release/functional_integration/Ssolidus_protein_function_master.tsv

QC_OUT=${RESULTS}/qc_clustering_no_mt_pca
DOUBLETS_OUT=${RESULTS}/doublet_reclustering_no_mt_pca
MARKERS_OUT=${RESULTS}/cluster_qc_markers_after_doublet_removal
VALIDATION_OUT=${RESULTS}/manual_cluster_validation
FINAL_OUT=${RESULTS}/final_cell_type_annotation

require_file() {
  [[ -s "$1" ]] || { echo "Missing file: $1" >&2; exit 1; }
}

require_dir() {
  [[ -d "$1" ]] || { echo "Missing directory: $1" >&2; exit 1; }
}

require_dir "${BIG_MATRIX}"
require_dir "${BREEDING_MATRIX}"
require_file "${FUNCTION_MASTER}"
mkdir -p "${QC_OUT}" "${DOUBLETS_OUT}" "${MARKERS_OUT}" "${VALIDATION_OUT}" "${FINAL_OUT}"

conda activate /work/cyu/conda_envs/ssol_scrna
unset R_LIBS R_LIBS_SITE
export R_LIBS_USER=${CONDA_PREFIX}/lib/R/library

Rscript --vanilla "${SCRIPTS}/03_seurat_qc_clustering_no_mt_pca.R" \
  "${BIG_MATRIX}" \
  "${BREEDING_MATRIX}" \
  "${QC_OUT}"

Rscript --vanilla "${SCRIPTS}/04_doublet_detection_reclustering.R" \
  "${QC_OUT}/ssolidus_nuclear_plus_mt_seurat.rds" \
  "${DOUBLETS_OUT}"

Rscript --vanilla "${SCRIPTS}/05_cluster_qc_and_markers.R" \
  "${DOUBLETS_OUT}/ssolidus_singlets_no_mt_pca_seurat.rds" \
  "${MARKERS_OUT}"

Rscript --vanilla "${SCRIPTS}/06_annotate_markers_with_genome_function.R" \
  "${MARKERS_OUT}/top20_markers_per_cluster.csv" \
  "${FUNCTION_MASTER}" \
  "${MARKERS_OUT}/top20_markers_per_cluster.functionally_annotated.tsv"

Rscript --vanilla "${SCRIPTS}/07_manual_cluster_validation.R" \
  "${MARKERS_OUT}/ssolidus_joined_layers_with_markers.rds" \
  "${FUNCTION_MASTER}" \
  "${VALIDATION_OUT}"

Rscript --vanilla "${SCRIPTS}/08_assign_final_cell_types.R" \
  "${MARKERS_OUT}/ssolidus_joined_layers_with_markers.rds" \
  "${FINAL_OUT}"

Rscript --vanilla "${SCRIPTS}/09_publication_figures.R" \
  "${FINAL_OUT}/ssolidus_final_annotated_seurat.rds" \
  "${MARKERS_OUT}/all_cluster_markers.csv" \
  "${FINAL_OUT}"

bash "${SCRIPTS}/10_record_analysis_manifest.sh"

echo "Final analysis completed"
echo "Annotated object: ${FINAL_OUT}/ssolidus_final_annotated_seurat.rds"
echo "Publication figures: ${FINAL_OUT}"
