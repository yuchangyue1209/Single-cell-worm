#!/usr/bin/env Rscript

# Calculate final cluster QC, sample composition and positive nuclear markers.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 05_cluster_qc_and_markers.R SINGLET_RDS OUTDIR")

input_rds <- args[1]
outdir <- args[2]
if (!file.exists(input_rds)) stop("Missing input RDS: ", input_rds)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
obj$seurat_clusters <- factor(as.character(obj$seurat_clusters), levels = as.character(0:24))
Idents(obj) <- "seurat_clusters"
DefaultAssay(obj) <- "RNA"
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

metadata <- obj[[]]
cluster_qc <- do.call(rbind, lapply(split(metadata, metadata$seurat_clusters), function(x) {
  data.frame(
    cluster = as.character(x$seurat_clusters[1]),
    cells = nrow(x),
    median_genes = median(x$nFeature_RNA),
    mean_genes = mean(x$nFeature_RNA),
    median_UMI = median(x$nCount_RNA),
    mean_UMI = mean(x$nCount_RNA),
    median_percent_mt = median(x$percent.mt),
    median_doublet_score = median(x$scDblFinder.score)
  )
}))
write.csv(cluster_qc, file.path(outdir, "cluster_QC_summary.csv"), row.names = FALSE)

composition <- as.data.frame(table(
  cluster = metadata$seurat_clusters,
  sample = metadata$sample
))
composition <- composition[composition$Freq > 0, , drop = FALSE]
composition$total_cluster_cells <- ave(composition$Freq, composition$cluster, FUN = sum)
composition$percent_within_cluster <- 100 * composition$Freq / composition$total_cluster_cells
write.csv(
  composition,
  file.path(outdir, "cluster_sample_composition.csv"),
  row.names = FALSE
)

message("Finding positive markers for all clusters")
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
markers <- markers[
  !grepl("^(nad[1-6]|nad4l|cox[1-3]|atp6|cob)$", markers$gene),
  ,
  drop = FALSE
]
write.csv(markers, file.path(outdir, "all_cluster_markers.csv"), row.names = FALSE)

top20 <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) {
  head(x[order(-x$avg_log2FC), , drop = FALSE], 20)
}))
write.csv(top20, file.path(outdir, "top20_markers_per_cluster.csv"), row.names = FALSE)

saveRDS(obj, file.path(outdir, "ssolidus_joined_layers_with_markers.rds"))

pdf(file.path(outdir, "umap_cluster_sample_split.pdf"), width = 12, height = 6)
print(DimPlot(obj, label = TRUE, repel = TRUE) + DimPlot(obj, split.by = "sample"))
dev.off()

message("Marker analysis completed")
message("Results: ", outdir)
