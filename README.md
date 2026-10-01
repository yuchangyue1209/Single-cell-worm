# *Schistocephalus solidus* single-cell RNA-seq atlas

This repository records the reproducible analysis of the Big and Breeding 10x libraries from Cell Ranger quantification through the final 25-cluster cell atlas. Commands are run on `stickleback`; `/work/cyu` paths are therefore not expected to exist on the desktop copy.

## Final analysis snapshot

- Reference: `/work/cyu/scRNA/reference_test/ssol_final_nuclear_plus_mt`
- Input runs: `Big_braker_tsa_utr_tagseq_plus_187_mt_2026` and `Breeding_braker_tsa_utr_tagseq_plus_187_mt_2026`
- QC-passing cells before doublet removal: 19,829
- Doublets: Big 658/8,766 (7.51%); Breeding 976/11,063 (8.82%)
- Final singlets: 18,195 (Big 8,108; Breeding 10,087)
- Final clustering: RNA SNN resolution 0.5, 25 clusters
- Top-marker rows matched to functional annotation: 481/500 (96.2%)
- Mitochondrial genes are retained for QC and expression analysis but excluded from variable-feature selection, scaling, PCA, neighbors, UMAP, and nuclear marker testing.

## Workflow

Run scripts in numerical order:

1. `01_build_nuclear_plus_mt_reference.sh`: combined Cell Ranger reference.
2. `02_cellranger_braker3.sh`: Big/Breeding quantification.
3. `03_seurat_qc_clustering_no_mt_pca.R`: sample QC and initial PCA.
4. `04_doublet_detection_reclustering.R`: per-sample scDblFinder and final clustering.
5. `05_cluster_qc_and_markers.R`: cluster QC and nuclear marker genes.
6. `06_annotate_markers_with_genome_function.R`: Swiss-Prot/eggNOG/InterPro integration.
7. `07_manual_cluster_validation.R`: targeted cluster comparisons, cluster-1 doublet-score audit, unresolved-cluster markers, and functional-module summaries used during manual review.
8. `08_assign_final_cell_types.R`: frozen detailed and broad labels for clusters 0–24.
9. `09_publication_figures.R`: final ordered UMAPs, descriptive sample composition, curated marker DotPlot, and numeric-order top-three-marker heatmap.
10. `10_record_analysis_manifest.sh`: file checksums, software versions, R session information, and the final numerical snapshot.

The exact full commands are in `run_final_pipeline.sh`. The script uses the real Cell Ranger run names, final functional-annotation table, and final output paths from this project. Existing output directories are never removed.

To regenerate only the downstream analysis from the completed 18,195-singlet object, use `run_from_existing_singlets.sh`. This avoids rerunning Cell Ranger and scDblFinder.

The Big and Breeding libraries are biological stage/sample groups without replicate libraries. Sample composition and stage-specific patterns are therefore descriptive; they are not formal replicated differential-abundance tests.

## Final outputs

```text
/work/cyu/scRNA/downstream/results/doublet_reclustering_no_mt_pca/
/work/cyu/scRNA/downstream/results/cluster_qc_markers_after_doublet_removal/
/work/cyu/scRNA/downstream/results/final_nuclear_markers_after_doublet_removal/
/work/cyu/scRNA/downstream/results/manual_cluster_validation/
/work/cyu/scRNA/downstream/results/final_cell_type_annotation/
```

Final annotated object:

```text
/work/cyu/scRNA/downstream/results/final_cell_type_annotation/ssolidus_final_annotated_seurat.rds
```

`scrna.sh` contains early STARsolo/reference experiments and is retained only for provenance.

## Recommended restart command

After copying this repository to `/work/cyu/scRNA/downstream`, run:

```bash
bash /work/cyu/scRNA/downstream/run_from_existing_singlets.sh
```

Use `run_final_pipeline.sh` only when the analysis must be rebuilt from the Cell Ranger matrices.

Removal of the 6,055-bp host-derived `scaffold_61` does not require Cell Ranger to be rerun: the scaffold contained no annotated genes and was already absent from the nuclear FASTA used in the final checks.
