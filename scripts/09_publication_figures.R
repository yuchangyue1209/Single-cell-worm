#!/usr/bin/env Rscript

# Generate the final publication figures used for the S. solidus cell atlas.
#
# Usage:
#   Rscript 09_publication_figures.R \
#     ssolidus_final_annotated_seurat.rds \
#     all_cluster_markers.csv \
#     output_directory

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: 09_publication_figures.R FINAL_RDS ALL_MARKERS_CSV OUTDIR",
    call. = FALSE
  )
}

rds_file <- args[1]
marker_file <- args[2]
outdir <- args[3]

if (!file.exists(rds_file)) stop("Missing RDS: ", rds_file)
if (!file.exists(marker_file)) stop("Missing marker table: ", marker_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading final annotated object")
obj <- readRDS(rds_file)
DefaultAssay(obj) <- "RNA"

if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

cluster_levels <- as.character(0:24)
obj$seurat_clusters <- factor(
  as.character(obj$seurat_clusters),
  levels = cluster_levels
)
Idents(obj) <- "seurat_clusters"

required_metadata <- c("sample", "cell_type", "broad_class")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0) {
  stop("Missing metadata: ", paste(missing_metadata, collapse = ", "))
}

save_plot <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(outdir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

# Figure 1A: numbered clusters in strict numeric order.
p1a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = FALSE
) +
  ggtitle("Cell clusters") +
  labs(color = "Cluster") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

save_plot("Fig1A_numbered_cluster_UMAP.pdf", p1a, 9, 7)

# Figure 1B: broad biological classes.
p1b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_class",
  label = TRUE,
  repel = TRUE,
  label.size = 3.2,
  raster = FALSE
) +
  ggtitle("Schistocephalus solidus cell atlas") +
  labs(color = "Broad cell class") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

save_plot("Fig1B_broad_cell_class_UMAP.pdf", p1b, 10, 7)

# Figure 1C: same embedding split by Big and Breeding samples.
p1c <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_class",
  split.by = "sample",
  ncol = 2,
  raster = FALSE
) +
  ggtitle("Broad cell classes by sample/stage") +
  labs(color = "Broad cell class") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_plot("Fig1C_broad_classes_by_sample.pdf", p1c, 14, 7)

# Figure 1D: descriptive cell composition. There are no biological replicates,
# so this is not used as an inferential differential-abundance test.
composition <- as.data.frame(table(
  sample = obj$sample,
  broad_class = obj$broad_class
))
composition <- composition[composition$Freq > 0, , drop = FALSE]
composition$total <- ave(composition$Freq, composition$sample, FUN = sum)
composition$percent <- 100 * composition$Freq / composition$total

write.csv(
  composition,
  file.path(outdir, "Fig1D_broad_class_composition.csv"),
  row.names = FALSE
)

p1d <- ggplot(
  composition,
  aes(x = sample, y = percent, fill = broad_class)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Broad cell-class composition",
    x = NULL,
    y = "Cells within sample (%)",
    fill = "Broad cell class"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

save_plot("Fig1D_broad_class_composition.pdf", p1d, 8, 6)

# Supplementary detailed annotation UMAPs.
p_s1a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 2.8,
  raster = FALSE
) +
  ggtitle("Schistocephalus solidus cell atlas") +
  labs(color = "Cell type") +
  theme_classic(base_size = 11)

save_plot("FigS1A_detailed_cell_type_UMAP.pdf", p_s1a, 14, 8)

p_s1b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "sample",
  ncol = 2,
  raster = FALSE
) +
  ggtitle("Detailed cell types by sample/stage") +
  labs(color = "Cell type") +
  theme_classic(base_size = 10)

save_plot("FigS1B_detailed_cell_types_by_sample.pdf", p_s1b, 16, 8)

# Figure 2A: manually curated marker panel used in the final interpretation.
# Gene IDs below reproduce the final panel documented in the analysis record.
marker_panel <- list(
  "Tegumental / epithelial" = c("g2683", "g611", "g3035", "g1496"),
  "Digestive" = c("g49", "g6117", "g9619"),
  "Metabolic / secretory" = c("g4444", "g8569", "g409"),
  "Muscle" = c("g6749", "g3376", "g2001"),
  "Neural" = c("g979", "g8158", "g7989", "g4206"),
  "Ciliated / sensory" = c("g4849", "g4850", "g5786"),
  "Protonephridial" = c("g1762", "g3", "g5045"),
  "Reproductive" = c("g7474", "g8058", "g1273", "g5195"),
  "Progenitor" = c("g1085", "g2132", "g1855"),
  "Rare neurosecretory" = c("g2648", "g8920", "g6919")
)

panel_table <- do.call(
  rbind,
  lapply(names(marker_panel), function(group) {
    data.frame(
      broad_class = group,
      gene = marker_panel[[group]],
      stringsAsFactors = FALSE
    )
  })
)

panel_table$present_in_object <- panel_table$gene %in% rownames(obj)
write.table(
  panel_table,
  file.path(outdir, "candidate_publication_marker_panel.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

features <- panel_table$gene[panel_table$present_in_object]
feature_groups <- panel_table$gene[panel_table$present_in_object]
names(feature_groups) <- panel_table$broad_class[panel_table$present_in_object]

if (length(features) > 0) {
  p2a <- DotPlot(
    obj,
    features = split(features, names(feature_groups)),
    group.by = "seurat_clusters",
    cols = c("grey90", "blue"),
    dot.scale = 7
  ) +
    RotatedAxis() +
    ggtitle("Selected markers supporting cell-cluster annotation") +
    labs(x = "Marker genes", y = "Cluster") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 8)
    )

  save_plot("Fig2A_representative_marker_DotPlot.pdf", p2a, 18, 10)
}

# Figure 2B: top three positive nuclear markers per cluster.
markers <- read.csv(marker_file, check.names = FALSE)
required_marker_columns <- c("cluster", "gene", "avg_log2FC")
if (!all(required_marker_columns %in% colnames(markers))) {
  stop("Marker table lacks: ", paste(setdiff(required_marker_columns, colnames(markers)), collapse = ", "))
}

markers$cluster <- factor(as.character(markers$cluster), levels = cluster_levels)
markers <- markers[
  !is.na(markers$cluster) &
    markers$gene %in% rownames(obj) &
    !grepl("^(nad[1-6]|nad4l|cox[1-3]|atp6|cob)$", markers$gene),
  ,
  drop = FALSE
]
markers <- markers[order(markers$cluster, -markers$avg_log2FC), , drop = FALSE]

top3 <- do.call(
  rbind,
  lapply(split(markers, markers$cluster, drop = TRUE), function(x) {
    head(x[!duplicated(x$gene), , drop = FALSE], 3)
  })
)

top3_genes <- unique(top3$gene)
write.csv(top3, file.path(outdir, "Fig2B_top3_markers.csv"), row.names = FALSE)

if (length(top3_genes) > 0) {
  obj <- ScaleData(obj, features = top3_genes, verbose = FALSE)
  p2b <- DoHeatmap(
    obj,
    features = top3_genes,
    group.by = "seurat_clusters",
    raster = TRUE,
    size = 3
  ) +
    ggtitle("Top nuclear markers for each cell cluster") +
    theme(plot.title = element_text(face = "bold"))

  save_plot("Fig2B_top3_marker_heatmap.pdf", p2b, 15, 10)
}

writeLines(
  c(
    paste("Input RDS:", normalizePath(rds_file)),
    paste("Marker file:", normalizePath(marker_file)),
    paste("Cells:", ncol(obj)),
    paste("Features:", nrow(obj)),
    paste("Clusters:", length(unique(obj$seurat_clusters))),
    paste("Missing detailed labels:", sum(is.na(obj$cell_type))),
    paste("Missing broad labels:", sum(is.na(obj$broad_class))),
    "",
    capture.output(sessionInfo())
  ),
  file.path(outdir, "publication_figure_provenance.txt")
)

message("Publication figures completed")
message("Output: ", outdir)
