#!/usr/bin/env Rscript

# Recreate the original two-panel UMAP layout with the final singlet atlas.
#
# Usage:
#   Rscript 12_final_umap_sample_and_cluster.R \
#     FINAL_ANNOTATED_RDS \
#     OUTPUT_DIRECTORY

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: 12_final_umap_sample_and_cluster.R FINAL_ANNOTATED_RDS OUTPUT_DIRECTORY",
    call. = FALSE
  )
}

rds_file <- args[1]
outdir <- args[2]

if (!file.exists(rds_file)) stop("Missing input RDS: ", rds_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading final annotated Seurat object")
obj <- readRDS(rds_file)

required_metadata <- c("sample", "seurat_clusters")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0) {
  stop("Missing metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!"umap" %in% names(obj@reductions)) stop("The object has no UMAP reduction")

# Preserve the intended sample and numerical cluster order.
obj$sample <- factor(as.character(obj$sample), levels = c("Big", "Breeding"))
obj$seurat_clusters <- factor(
  as.character(obj$seurat_clusters),
  levels = as.character(0:24)
)
Idents(obj) <- "seurat_clusters"

if (anyNA(obj$sample)) stop("Unexpected sample label found")
if (anyNA(obj$seurat_clusters)) stop("Unexpected cluster label found")

sample_colors <- c("Big" = "#F8766D", "Breeding" = "#00BFC4")

p_sample <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "sample",
  cols = sample_colors,
  pt.size = 0.25,
  shuffle = TRUE,
  raster = FALSE
) +
  ggtitle("UMAP by sample/stage") +
  labs(color = NULL) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

p_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4.2,
  pt.size = 0.25,
  raster = FALSE
) +
  ggtitle("UMAP by cluster") +
  labs(color = NULL) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

combined <- p_sample + p_cluster +
  plot_layout(widths = c(1, 1), guides = "keep")

pdf_file <- file.path(outdir, "final_UMAP_by_sample_and_cluster_old_style.pdf")
png_file <- file.path(outdir, "final_UMAP_by_sample_and_cluster_old_style.png")

ggsave(
  pdf_file,
  combined,
  width = 16,
  height = 7.2,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  png_file,
  combined,
  width = 16,
  height = 7.2,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

cat("Cells:", ncol(obj), "\n")
cat("Clusters:", nlevels(droplevels(obj$seurat_clusters)), "\n")
cat("PDF:", pdf_file, "\n")
cat("PNG:", png_file, "\n")
