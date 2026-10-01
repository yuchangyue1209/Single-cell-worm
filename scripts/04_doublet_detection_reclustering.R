#!/usr/bin/env Rscript

# Detect doublets separately in each sample using scDblFinder, retain singlets,
# and reconstruct the final no-mitochondrial-feature PCA/UMAP graph.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 04_doublet_detection_reclustering.R INPUT_RDS OUTDIR")
}

input_rds <- args[1]
outdir <- args[2]
if (!file.exists(input_rds)) stop("Missing input RDS: ", input_rds)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

obj <- readRDS(input_rds)
sample_objects <- SplitObject(obj, split.by = "sample")

run_scdblfinder <- function(sample_object) {
  sce <- as.SingleCellExperiment(sample_object)
  sce <- scDblFinder(sce)
  sample_object$scDblFinder.score <- colData(sce)$scDblFinder.score
  sample_object$scDblFinder.class <- colData(sce)$scDblFinder.class
  sample_object
}

message("Running scDblFinder independently for each sample")
sample_objects <- lapply(sample_objects, run_scdblfinder)
classified <- merge(sample_objects[[1]], y = sample_objects[-1])

doublet_summary <- as.data.frame(table(
  sample = classified$sample,
  classification = classified$scDblFinder.class
))
doublet_summary$total_cells <- ave(doublet_summary$Freq, doublet_summary$sample, FUN = sum)
doublet_summary$percent <- 100 * doublet_summary$Freq / doublet_summary$total_cells
write.csv(
  doublet_summary,
  file.path(outdir, "doublet_summary_by_sample.csv"),
  row.names = FALSE
)

message("Retaining singlets and rebuilding the dimensional reduction")
singlets <- subset(classified, subset = scDblFinder.class == "singlet")
singlets <- NormalizeData(singlets, verbose = FALSE)
singlets <- FindVariableFeatures(singlets, nfeatures = 2000, verbose = FALSE)
mt_genes <- grep(
  "^(nad[1-6]|nad4l|cox[1-3]|atp6|cob)$",
  rownames(singlets),
  value = TRUE
)
VariableFeatures(singlets) <- setdiff(VariableFeatures(singlets), mt_genes)
singlets <- ScaleData(singlets, features = VariableFeatures(singlets), verbose = FALSE)
singlets <- RunPCA(
  singlets,
  features = VariableFeatures(singlets),
  npcs = 50,
  verbose = FALSE
)
singlets <- FindNeighbors(singlets, dims = 1:30, verbose = FALSE)

resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0)
singlets <- FindClusters(singlets, resolution = resolutions, verbose = FALSE)
singlets$seurat_clusters <- singlets$RNA_snn_res.0.5
Idents(singlets) <- "seurat_clusters"
singlets <- RunUMAP(singlets, dims = 1:30, seed.use = 2026, verbose = FALSE)

resolution_columns <- paste0("RNA_snn_res.", resolutions)
resolution_counts <- data.frame(
  resolution_column = resolution_columns,
  clusters = vapply(resolution_columns, function(column) {
    length(unique(singlets[[column]][, 1]))
  }, integer(1))
)
write.csv(
  resolution_counts,
  file.path(outdir, "resolution_cluster_counts.csv"),
  row.names = FALSE
)

saveRDS(
  singlets,
  file.path(outdir, "ssolidus_singlets_no_mt_pca_seurat.rds")
)

pdf(file.path(outdir, "singlet_UMAP_clusters.pdf"), width = 8, height = 7)
print(DimPlot(singlets, label = TRUE, repel = TRUE))
dev.off()

message("Doublet removal and reclustering completed")
message("Singlet cells: ", ncol(singlets))
