#!/usr/bin/env Rscript

# Freeze the final detailed and broad cell-type labels for clusters 0-24.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 08_assign_final_cell_types.R INPUT_RDS OUTDIR")

input_rds <- args[1]
outdir <- args[2]
if (!file.exists(input_rds)) stop("Missing input RDS: ", input_rds)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

detailed_labels <- c(
  "0" = "Absorptive/tegumental cells",
  "1" = "Calcium-binding metabolic cells",
  "2" = "Low-RNA transitional cells",
  "3" = "Secretory/tegumental cells",
  "4" = "Proliferative vitellogenic cells",
  "5" = "Low-RNA stressed cells",
  "6" = "Contractile muscle",
  "7" = "Differentiating spermatogenic cells",
  "8" = "Proteolytic digestive cells",
  "9" = "Cycling germline/spermatogenic precursors",
  "10" = "Ciliated absorptive/digestive cells",
  "11" = "Neural/ECM-associated cells",
  "12" = "Glutamatergic neurosecretory cells",
  "13" = "Metabolic secretory/absorptive cells",
  "14" = "FoxA/GATA-positive epithelial cells",
  "15" = "ECM-associated contractile muscle",
  "16" = "Neural cells",
  "17" = "Neural cells II",
  "18" = "Flame/protonephridial cells",
  "19" = "Digestive/tegumental secretory cells",
  "20" = "Proteolytic digestive cells II",
  "21" = "Ciliated sensory/tegumental cells",
  "22" = "Cycling somatic/secretory progenitors",
  "23" = "Cycling secretory/tegumental progenitors",
  "24" = "Protease-inhibitor neurosecretory cells"
)

broad_labels <- c(
  "0" = "Tegumental/epithelial", "1" = "Metabolic/secretory",
  "2" = "Low-RNA/uncertain", "3" = "Tegumental/epithelial",
  "4" = "Reproductive", "5" = "Low-RNA/uncertain",
  "6" = "Muscle", "7" = "Reproductive", "8" = "Digestive",
  "9" = "Reproductive", "10" = "Digestive", "11" = "Neural",
  "12" = "Neural", "13" = "Metabolic/secretory",
  "14" = "Tegumental/epithelial", "15" = "Muscle",
  "16" = "Neural", "17" = "Neural", "18" = "Protonephridial",
  "19" = "Digestive", "20" = "Digestive",
  "21" = "Ciliated/sensory", "22" = "Progenitor",
  "23" = "Progenitor", "24" = "Neural"
)

clusters <- as.character(obj$seurat_clusters)
unknown <- setdiff(unique(clusters), names(detailed_labels))
if (length(unknown) > 0) stop("Unknown cluster IDs: ", paste(unknown, collapse = ", "))

obj$cell_type <- factor(
  unname(detailed_labels[clusters]),
  levels = unname(detailed_labels)
)
obj$broad_class <- factor(
  unname(broad_labels[clusters]),
  levels = unique(unname(broad_labels))
)
obj$seurat_clusters <- factor(clusters, levels = as.character(0:24))
Idents(obj) <- "seurat_clusters"

if (anyNA(obj$cell_type) || anyNA(obj$broad_class)) {
  stop("Some cells did not receive final labels")
}

saveRDS(obj, file.path(outdir, "ssolidus_final_annotated_seurat.rds"))

cell_counts <- as.data.frame(table(
  cluster = obj$seurat_clusters,
  cell_type = obj$cell_type
))
cell_counts <- cell_counts[cell_counts$Freq > 0, , drop = FALSE]
colnames(cell_counts)[3] <- "cells"
write.csv(cell_counts, file.path(outdir, "final_cell_type_counts.csv"), row.names = FALSE)

composition <- as.data.frame(table(
  sample = obj$sample,
  broad_class = obj$broad_class
))
composition <- composition[composition$Freq > 0, , drop = FALSE]
composition$total <- ave(composition$Freq, composition$sample, FUN = sum)
composition$percent <- 100 * composition$Freq / composition$total
write.csv(
  composition,
  file.path(outdir, "broad_cell_class_composition.csv"),
  row.names = FALSE
)

pdf(file.path(outdir, "final_cell_type_UMAP.pdf"), width = 13, height = 8)
print(DimPlot(
  obj, group.by = "cell_type", label = TRUE, repel = TRUE,
  label.size = 3, raster = FALSE
))
dev.off()

pdf(file.path(outdir, "final_cell_type_UMAP_by_sample.pdf"), width = 15, height = 7)
print(DimPlot(obj, group.by = "broad_class", split.by = "sample", raster = FALSE))
dev.off()

cat("Cells:", ncol(obj), "\n")
cat("Clusters:", length(unique(obj$seurat_clusters)), "\n")
cat("Missing labels:", sum(is.na(obj$cell_type)), "\n")
cat("Results:", outdir, "\n")
