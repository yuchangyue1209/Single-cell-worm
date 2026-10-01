#!/usr/bin/env bash
set -euo pipefail

# Record checksums, versions and the final analysis snapshot.

PROJECT=/work/cyu/scRNA
RESULTS=${PROJECT}/downstream/results
OUT=${RESULTS}/final_cell_type_annotation/reproducibility_manifest
FINAL_RDS=${RESULTS}/final_cell_type_annotation/ssolidus_final_annotated_seurat.rds
MARKERS=${RESULTS}/cluster_qc_markers_after_doublet_removal/all_cluster_markers.csv
CELL_TYPES=${RESULTS}/final_cell_type_annotation/final_cell_type_counts.csv
COMPOSITION=${RESULTS}/final_cell_type_annotation/Fig1D_broad_class_composition.csv
FUNCTION_MASTER=/work/cyu/annotation/ssol-annotation/results/final_filtered_release/functional_integration/Ssolidus_protein_function_master.tsv

mkdir -p "${OUT}"
for file in "${FINAL_RDS}" "${MARKERS}" "${CELL_TYPES}" "${FUNCTION_MASTER}"; do
  [[ -s "${file}" ]] || { echo "Missing file: ${file}" >&2; exit 1; }
done

sha256sum "${FINAL_RDS}" "${MARKERS}" "${CELL_TYPES}" "${FUNCTION_MASTER}" \
  > "${OUT}/SHA256SUMS"
if [[ -s "${COMPOSITION}" ]]; then
  sha256sum "${COMPOSITION}" >> "${OUT}/SHA256SUMS"
fi

{
  date -u '+generated_utc\t%Y-%m-%dT%H:%M:%SZ'
  printf 'host\t%s\n' "$(hostname)"
  printf 'conda_environment\t%s\n' "${CONDA_PREFIX:-not_recorded}"
  printf 'R\t'
  R --version | head -1
  printf 'cellranger\t'
  cellranger --version 2>/dev/null || echo 'not_on_current_PATH'
} > "${OUT}/software_versions.tsv"

Rscript --vanilla - "${OUT}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
packages <- c(
  "Seurat", "SeuratObject", "SingleCellExperiment", "scDblFinder",
  "Matrix", "ggplot2", "patchwork"
)
versions <- data.frame(
  package = packages,
  version = vapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(packageVersion(pkg))
    } else {
      "not_installed"
    }
  }, character(1))
)
write.table(
  versions, file.path(outdir, "R_package_versions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
RS

cat > "${OUT}/analysis_snapshot.tsv" <<'EOF'
metric	value
QC_passing_cells_before_doublet_removal	19829
Big_doublets	658
Big_total_before_doublet_removal	8766
Breeding_doublets	976
Breeding_total_before_doublet_removal	11063
final_singlets	18195
Big_final_singlets	8108
Breeding_final_singlets	10087
final_clusters	25
cluster_resolution	0.5
nuclear_genes_in_final_object	9498
marker_function_match	481/500 (96.2%)
EOF

echo "Manifest written to: ${OUT}"
