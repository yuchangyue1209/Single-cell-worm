#!/usr/bin/env Rscript

# Export four separate, publication-ready Figure 1 panels from the final
# 18,195-cell S. solidus atlas.
#
# Usage:
#   Rscript 13_nature_figure1_panels.R FINAL_ANNOTATED_RDS OUTPUT_DIRECTORY

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: 13_nature_figure1_panels.R FINAL_ANNOTATED_RDS OUTPUT_DIRECTORY",
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

cluster_levels <- as.character(0:24)
broad_levels <- c(
  "Tegumental/epithelial",
  "Digestive",
  "Metabolic/secretory",
  "Muscle",
  "Neural",
  "Ciliated/sensory",
  "Protonephridial",
  "Reproductive",
  "Progenitor",
  "Low-RNA/uncertain"
)

# Some saved versions of the final object contain the final cluster labels but
# predate creation of the broad_class metadata column. Reconstruct it from the
# frozen, manually validated cluster-to-class mapping when necessary.
broad_class_by_cluster <- c(
  "0" = "Tegumental/epithelial",
  "1" = "Metabolic/secretory",
  "2" = "Low-RNA/uncertain",
  "3" = "Tegumental/epithelial",
  "4" = "Reproductive",
  "5" = "Low-RNA/uncertain",
  "6" = "Muscle",
  "7" = "Reproductive",
  "8" = "Digestive",
  "9" = "Reproductive",
  "10" = "Digestive",
  "11" = "Neural",
  "12" = "Neural",
  "13" = "Metabolic/secretory",
  "14" = "Tegumental/epithelial",
  "15" = "Muscle",
  "16" = "Neural",
  "17" = "Neural",
  "18" = "Protonephridial",
  "19" = "Digestive",
  "20" = "Digestive",
  "21" = "Ciliated/sensory",
  "22" = "Progenitor",
  "23" = "Progenitor",
  "24" = "Neural"
)

cluster_values <- as.character(obj$seurat_clusters)
if (!"broad_class" %in% colnames(obj[[]])) {
  message("Metadata column broad_class is absent; reconstructing it from final clusters")
  unknown_clusters <- setdiff(unique(cluster_values), names(broad_class_by_cluster))
  if (length(unknown_clusters) > 0) {
    stop(
      "Cannot reconstruct broad_class for cluster(s): ",
      paste(unknown_clusters, collapse = ", ")
    )
  }
  obj$broad_class <- unname(broad_class_by_cluster[cluster_values])
}

obj$sample <- factor(as.character(obj$sample), levels = c("Big", "Breeding"))
obj$seurat_clusters <- factor(
  cluster_values,
  levels = cluster_levels
)
obj$broad_class <- factor(as.character(obj$broad_class), levels = broad_levels)
Idents(obj) <- "seurat_clusters"

if (anyNA(obj$sample)) stop("Unexpected sample label found")
if (anyNA(obj$seurat_clusters)) stop("Unexpected cluster label found")
if (anyNA(obj$broad_class)) stop("Unexpected broad_class label found")

# Fixed, colorblind-conscious palette used consistently in panels b-d.
broad_colors <- c(
  "Tegumental/epithelial" = "#D55E00",
  "Digestive" = "#E69F00",
  "Metabolic/secretory" = "#A6761D",
  "Muscle" = "#009E73",
  "Neural" = "#0072B2",
  "Ciliated/sensory" = "#56B4E9",
  "Protonephridial" = "#CC79A7",
  "Reproductive" = "#C51B7D",
  "Progenitor" = "#6A3D9A",
  "Low-RNA/uncertain" = "#8C8C8C"
)

# A reproducible 25-color qualitative palette for numerical clusters.
cluster_colors <- setNames(hcl.colors(25, palette = "Dynamic"), cluster_levels)

# Nature recommends a consistent sans-serif typeface and 5-7 pt lettering at
# final size; panel letters are 8 pt bold. The panels below are saved close to
# their intended printed dimensions so the point sizes are meaningful.
nature_theme <- function(panel_letter) {
  theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.title = element_text(size = 7, color = "black"),
      axis.text = element_text(size = 6, color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold", color = "black"),
      plot.tag.position = c(0, 1),
      legend.title = element_text(size = 6.5, face = "bold"),
      legend.text = element_text(size = 5.5),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(3.2, "mm"),
      legend.spacing.y = unit(0.2, "mm"),
      plot.margin = margin(5, 5, 4, 5, unit = "mm")
    ) +
    labs(tag = panel_letter, x = "UMAP 1", y = "UMAP 2")
}

save_panel <- function(plot, stem, width_mm, height_mm) {
  pdf_file <- file.path(outdir, paste0(stem, ".pdf"))
  tif_file <- file.path(outdir, paste0(stem, ".tiff"))
  png_file <- file.path(outdir, paste0(stem, ".png"))

  ggsave(
    pdf_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, limitsize = FALSE
  )
  ggsave(
    tif_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = 600, compression = "lzw", limitsize = FALSE
  )
  ggsave(
    png_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = 600, limitsize = FALSE
  )
}

# Seurat does not expose point transparency directly through DimPlot().
# Apply alpha only to point layers, leaving text labels fully opaque.
set_point_alpha <- function(plot, alpha = 0.70) {
  for (i in seq_along(plot$layers)) {
    geom_class <- class(plot$layers[[i]]$geom)
    if (any(geom_class %in% c("GeomPoint", "GeomScattermore"))) {
      plot$layers[[i]]$aes_params$alpha <- alpha
    }
  }
  plot
}

# Panel a: final numerical clusters in strict order 0-24.
p_a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = cluster_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 2.4,
  pt.size = 0.05,
  raster = FALSE
) +
  ggtitle("Cell clusters") +
  nature_theme("a") +
  theme(legend.position = "none")
p_a <- set_point_alpha(p_a, 0.70)

# Panel b: broad biological classes using the same palette as panels c-d.
p_b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_class",
  cols = broad_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 2.35,
  pt.size = 0.05,
  raster = FALSE
) +
  ggtitle("Broad cell classes") +
  nature_theme("b") +
  theme(legend.position = "none")
p_b <- set_point_alpha(p_b, 0.70)

# Panel c: broad classes displayed separately for each sample/stage.
p_c <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_class",
  split.by = "sample",
  cols = adjustcolor(broad_colors, alpha.f = 0.70),
  ncol = 2,
  pt.size = 0.045,
  raster = FALSE
) +
  ggtitle("Broad cell classes by sample") +
  nature_theme("c") +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
    strip.text = element_text(size = 6.5, face = "plain")
  )
p_c <- p_c & theme(legend.position = "none")

# Panel d: descriptive composition only; no replicate-based inference is made.
composition <- as.data.frame(table(
  sample = obj$sample,
  broad_class = obj$broad_class
))
composition <- composition[composition$Freq > 0, , drop = FALSE]
composition$total <- ave(composition$Freq, composition$sample, FUN = sum)
composition$percent <- 100 * composition$Freq / composition$total

write.csv(
  composition,
  file.path(outdir, "Nature_Fig1D_broad_class_composition.csv"),
  row.names = FALSE
)

p_d <- ggplot(
  composition,
  aes(x = sample, y = percent, fill = broad_class)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.20) +
  scale_fill_manual(values = broad_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  ggtitle("Broad cell-class composition") +
  labs(x = NULL, y = "Cells within sample (%)", fill = "Broad cell class") +
  guides(fill = guide_legend(ncol = 1)) +
  nature_theme("d") +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

# Individual panels. Panels a, b and d are single-column width (89 mm).
# Panel c is wider because it contains two sample facets.
save_panel(p_a, "Nature_Fig1A_cluster_UMAP", 89, 82)
save_panel(p_b, "Nature_Fig1B_broad_cell_class_UMAP", 89, 82)
save_panel(p_c, "Nature_Fig1C_broad_classes_by_sample", 136, 82)
save_panel(p_d, "Nature_Fig1D_broad_class_composition", 89, 82)

cat("Cells:", ncol(obj), "\n")
cat("Clusters:", nlevels(droplevels(obj$seurat_clusters)), "\n")
cat("Output directory:", outdir, "\n")
cat("Generated four PDF, TIFF and PNG panels plus the panel-d source table.\n")
