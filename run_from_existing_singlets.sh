#!/usr/bin/env bash
set -euo pipefail

# Recommended restart point: regenerate downstream analyses from the completed
# 18,195-singlet Seurat object without rerunning Cell Ranger or scDblFinder.

PROJECT=/work/cyu/scRNA
SCRIPTS=${PROJECT}/downstream/scripts
RESULTS=${PROJECT}/downstream/results
FUNCTION_MASTER=/work/cyu/annotation/ssol-annotation/results/final_filtered_release/functional_integration/Ssolidus_protein_function_master.tsv
SINGLET_RDS=${RESULTS}/doublet_reclustering_no_mt_pca/ssolidus_singlets_no_mt_pca_seurat.rds
MARKERS_OUT=${RESULTS}/cluster_qc_markers_after_doublet_removal
VALIDATION_OUT=${RESULTS}/manual_cluster_validation
FINAL_OUT=${RESULTS}/final_cell_type_annotation

for file in "${SINGLET_RDS}" "${FUNCTION_MASTER}"; do
  [[ -s "${file}" ]] || { echo "Missing file: ${file}" >&2; exit 1; }
done
mkdir -p "${MARKERS_OUT}" "${VALIDATION_OUT}" "${FINAL_OUT}"

conda activate /work/cyu/conda_envs/ssol_scrna
unset R_LIBS R_LIBS_SITE
export R_LIBS_USER=${CONDA_PREFIX}/lib/R/library

Rscript --vanilla "${SCRIPTS}/05_cluster_qc_and_markers.R" \
  "${SINGLET_RDS}" "${MARKERS_OUT}"

Rscript --vanilla "${SCRIPTS}/06_annotate_markers_with_genome_function.R" \
  "${MARKERS_OUT}/top20_markers_per_cluster.csv" \
  "${FUNCTION_MASTER}" \
  "${MARKERS_OUT}/top20_markers_per_cluster.functionally_annotated.tsv"

Rscript --vanilla "${SCRIPTS}/07_manual_cluster_validation.R" \
  "${MARKERS_OUT}/ssolidus_joined_layers_with_markers.rds" \
  "${FUNCTION_MASTER}" "${VALIDATION_OUT}"

Rscript --vanilla "${SCRIPTS}/08_assign_final_cell_types.R" \
  "${MARKERS_OUT}/ssolidus_joined_layers_with_markers.rds" "${FINAL_OUT}"

Rscript --vanilla "${SCRIPTS}/09_publication_figures.R" \
  "${FINAL_OUT}/ssolidus_final_annotated_seurat.rds" \
  "${MARKERS_OUT}/all_cluster_markers.csv" "${FINAL_OUT}"

bash "${SCRIPTS}/10_record_analysis_manifest.sh"

echo "Analysis regenerated from the existing singlet object"
echo "Final results: ${FINAL_OUT}"
