#!/usr/bin/env Rscript

# Initial Seurat QC, normalization and clustering for Big and Breeding samples.
# Mitochondrial genes are retained in the expression matrix and used for QC,
# but excluded from variable features, scaling and PCA.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage: 03_seurat_qc_clustering_no_mt_pca.R",
    "BIG_MATRIX BREEDING_MATRIX OUTDIR",
    "[MIN_FEATURES=200] [MAX_FEATURES=6000] [MAX_MT=25]"
  ))
}

big_dir <- args[1]
breeding_dir <- args[2]
outdir <- args[3]
min_features <- if (length(args) >= 4) as.numeric(args[4]) else 200
max_features <- if (length(args) >= 5) as.numeric(args[5]) else 6000
max_mt <- if (length(args) >= 6) as.numeric(args[6]) else 25
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mitochondrial_pattern <- "^(nad[1-6]|nad4l|cox[1-3]|atp6|cob)$"

read_sample <- function(matrix_dir, sample_name) {
  if (!dir.exists(matrix_dir)) stop("Missing matrix directory: ", matrix_dir)
  counts <- Read10X(matrix_dir)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3
  )
  object$sample <- sample_name
  object[["percent.mt"]] <- PercentageFeatureSet(
    object,
    pattern = mitochondrial_pattern
  )
  object
}

message("Reading Cell Ranger matrices")
big <- read_sample(big_dir, "Big")
breeding <- read_sample(breeding_dir, "Breeding")

before <- data.frame(
  sample = c("Big", "Breeding"),
  cells = c(ncol(big), ncol(breeding))
)

filter_sample <- function(object) {
  subset(
    object,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      percent.mt <= max_mt
  )
}

big <- filter_sample(big)
breeding <- filter_sample(breeding)

after <- data.frame(
  sample = c("Big", "Breeding"),
  cells = c(ncol(big), ncol(breeding))
)
qc_counts <- merge(before, after, by = "sample", suffixes = c("_before", "_after"))
write.csv(qc_counts, file.path(outdir, "QC_cell_counts.csv"), row.names = FALSE)

message("Merging, normalizing and clustering")
obj <- merge(big, y = breeding, add.cell.ids = c("Big", "Breeding"))
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
mt_genes <- grep(mitochondrial_pattern, rownames(obj), value = TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj), mt_genes)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 2026, verbose = FALSE)

saveRDS(obj, file.path(outdir, "ssolidus_nuclear_plus_mt_seurat.rds"))
write.csv(obj[[]], file.path(outdir, "cell_metadata_after_QC.csv"))

pdf(file.path(outdir, "qc_violin_by_sample.pdf"), width = 10, height = 6)
print(VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
))
dev.off()

pdf(file.path(outdir, "umap_by_sample_and_cluster.pdf"), width = 12, height = 5)
print(DimPlot(obj, group.by = "sample") + DimPlot(obj, label = TRUE))
dev.off()

writeLines(c(
  paste("min_features", min_features),
  paste("max_features", max_features),
  paste("max_mt", max_mt),
  paste("mt_features_excluded_from_PCA", length(mt_genes)),
  "",
  capture.output(sessionInfo())
), file.path(outdir, "analysis_provenance.txt"))

message("Initial QC and clustering completed")
message("Results: ", outdir)
