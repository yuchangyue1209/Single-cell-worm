#!/usr/bin/env Rscript

# Reproduce the targeted cluster-validation analyses used to assign the final
# S. solidus cell atlas labels.
# Usage: Rscript 07_manual_cluster_validation.R INPUT_RDS FUNCTION_MASTER_TSV OUTDIR

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 07_manual_cluster_validation.R INPUT_RDS FUNCTION_MASTER_TSV OUTDIR")
}

rds_file <- args[1]
master_file <- args[2]
outdir <- args[3]
if (!file.exists(rds_file)) stop("Missing RDS: ", rds_file)
if (!file.exists(master_file)) stop("Missing functional table: ", master_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading Seurat object")
obj <- readRDS(rds_file)
DefaultAssay(obj) <- "RNA"
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
cluster_levels <- as.character(0:24)
obj$seurat_clusters <- factor(as.character(obj$seurat_clusters), levels = cluster_levels)
Idents(obj) <- "seurat_clusters"

message("Reading functional annotation")
functions <- read.delim(
  master_file, check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE
)
if (!"gene_id" %in% colnames(functions)) stop("Functional table lacks gene_id")
if ("evidence_source_count" %in% colnames(functions)) {
  functions <- functions[order(functions$gene_id, -functions$evidence_source_count), ]
} else {
  functions <- functions[order(functions$gene_id), ]
}
functions <- functions[!duplicated(functions$gene_id), , drop = FALSE]

annotation_columns <- intersect(
  c("product", "best_annotation", "swissprot_product", "eggnog_preferred_name",
    "eggnog_description", "interpro_descriptions"),
  colnames(functions)
)
if (length(annotation_columns) == 0) stop("No usable annotation columns")
functions$annotation_text <- apply(
  functions[, annotation_columns, drop = FALSE], 1,
  function(x) paste(unique(x[!is.na(x) & nzchar(x) & x != "-"]), collapse = "; ")
)

write_tsv <- function(x, filename) {
  write.table(x, file.path(outdir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
}

annotate_markers <- function(markers) {
  markers$gene <- rownames(markers)
  markers$pct_difference <- markers$pct.1 - markers$pct.2
  markers$.row_order <- seq_len(nrow(markers))
  out <- merge(markers, functions, by.x = "gene", by.y = "gene_id",
               all.x = TRUE, sort = FALSE)
  out <- out[order(out$.row_order), , drop = FALSE]
  out$.row_order <- NULL
  out
}

find_positive_markers <- function(cluster, min_pct = 0.05, logfc = 0.25) {
  markers <- FindMarkers(
    obj, ident.1 = as.character(cluster), only.pos = TRUE,
    min.pct = min_pct, logfc.threshold = logfc, test.use = "wilcox"
  )
  markers <- markers[order(-markers$avg_log2FC), , drop = FALSE]
  annotate_markers(markers)
}

write_cluster_qc <- function(clusters, prefix) {
  metadata <- obj[[]]
  metadata$cluster <- as.character(obj$seurat_clusters)
  selected <- metadata[metadata$cluster %in% as.character(clusters), , drop = FALSE]
  composition <- as.data.frame(table(cluster = selected$cluster, sample = selected$sample))
  composition <- composition[composition$Freq > 0, , drop = FALSE]
  composition$total <- ave(composition$Freq, composition$cluster, FUN = sum)
  composition$percent <- 100 * composition$Freq / composition$total
  write.csv(composition, file.path(outdir, paste0(prefix, "_sample_composition.csv")),
            row.names = FALSE)

  qc <- do.call(rbind, lapply(split(selected, selected$cluster), function(x) {
    data.frame(
      cluster = x$cluster[1], cells = nrow(x),
      median_genes = median(x$nFeature_RNA), median_UMI = median(x$nCount_RNA),
      median_percent_mt = median(x$percent.mt),
      median_doublet_score = if ("scDblFinder.score" %in% colnames(x)) {
        median(x$scDblFinder.score)
      } else NA_real_
    )
  }))
  write.csv(qc, file.path(outdir, paste0(prefix, "_QC_summary.csv")), row.names = FALSE)
}

# Priority comparisons used in the final interpretation.
write_cluster_qc(c(6, 15), "cluster6_cluster15")
cluster15_vs_6 <- FindMarkers(
  obj, ident.1 = "15", ident.2 = "6", only.pos = FALSE,
  min.pct = 0.05, logfc.threshold = 0.25, test.use = "wilcox"
)
cluster15_vs_6 <- annotate_markers(cluster15_vs_6)
cluster15_vs_6 <- cluster15_vs_6[order(-cluster15_vs_6$avg_log2FC), ]
write_tsv(cluster15_vs_6[cluster15_vs_6$avg_log2FC > 0, , drop = FALSE],
          "cluster15_top_specific_markers.tsv")
write_tsv(cluster15_vs_6[cluster15_vs_6$avg_log2FC < 0, , drop = FALSE],
          "cluster6_top_specific_markers.tsv")

write_cluster_qc(c(7, 9), "cluster7_cluster9")
cluster7_vs_9 <- FindMarkers(
  obj, ident.1 = "7", ident.2 = "9", only.pos = FALSE,
  min.pct = 0.05, logfc.threshold = 0.25, test.use = "wilcox"
)
cluster7_vs_9 <- annotate_markers(cluster7_vs_9)
cluster7_vs_9 <- cluster7_vs_9[order(-cluster7_vs_9$avg_log2FC), ]
write_tsv(cluster7_vs_9, "cluster7_vs_cluster9_markers.annotated.tsv")

write_cluster_qc(4, "cluster4")
cluster4 <- find_positive_markers(4)
write_tsv(cluster4, "cluster4_strong_markers.annotated.tsv")
oogenic_pattern <- paste(c(
  "oocyt", "oogen", "vitell", "egg antigen", "germ cell", "meiotic",
  "meiosis", "CPEB", "polyadenylation element", "maternal", "embryo"
), collapse = "|")
write_tsv(
  cluster4[grepl(oogenic_pattern, cluster4$annotation_text, ignore.case = TRUE), ],
  "cluster4_oogenic_keyword_hits.tsv"
)

write_cluster_qc(c(22, 23), "cluster22_cluster23")
for (cluster in c(22, 23)) {
  write_tsv(find_positive_markers(cluster),
            paste0("cluster", cluster, "_progenitor_markers.annotated.tsv"))
}

write_cluster_qc(24, "cluster24")
write_tsv(find_positive_markers(24, 0.01, 0.1), "cluster24_markers.annotated.tsv")

# Cluster 1 doublet-score audit.
if ("scDblFinder.score" %in% colnames(obj[[]])) {
  cells1 <- WhichCells(obj, idents = "1")
  cutoff <- unname(quantile(obj$scDblFinder.score[cells1], 0.75, na.rm = TRUE))
  obj1 <- subset(obj, cells = cells1)
  obj1$cluster1_score_group <- ifelse(
    obj1$scDblFinder.score >= cutoff, "upper_25pct_score", "lower_75pct_score"
  )
  Idents(obj1) <- "cluster1_score_group"
  score_summary <- aggregate(
    cbind(nFeature_RNA, nCount_RNA, percent.mt, scDblFinder.score) ~ cluster1_score_group,
    data = obj1[[]], FUN = median
  )
  write.csv(score_summary, file.path(outdir, "cluster1_doublet_score_group_QC.csv"),
            row.names = FALSE)
  score_markers <- FindMarkers(
    obj1, ident.1 = "upper_25pct_score", ident.2 = "lower_75pct_score",
    only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  write_tsv(annotate_markers(score_markers),
            "cluster1_upper25_vs_lower75_doublet_score_markers.tsv")
}

# Full marker review for the previously unresolved clusters.
unresolved_clusters <- c(0, 1, 2, 3, 5, 10, 13, 14, 21)
write_cluster_qc(unresolved_clusters, "unresolved_cluster")
for (cluster in unresolved_clusters) {
  message("Processing unresolved cluster ", cluster)
  write_tsv(head(find_positive_markers(cluster), 50),
            paste0("cluster", cluster, "_top50_markers.annotated.tsv"))
}

# Functional module scores are supporting evidence, not stand-alone labels.
program_patterns <- list(
  tegument_membrane = "tegument|tetraspanin|epithelial|membrane|surface protein",
  secretory_proteolytic = "secret|protease|peptidase|cathepsin|trypsin|chymotrypsin",
  digestive_epithelial = "digest|absorptive|lysosom|aminopeptidase|lipase",
  calcium_contractile = "myosin|actin|troponin|tropomyosin|calponin|contract|calcium-binding",
  neural_sensory = "neural|neuron|synap|neuropeptide|glutamat|acetylcholine|sensory",
  cilia_microtubule = "cilia|ciliary|dynein|intraflagellar|microtubule|tubulin",
  reproductive_oogenic = "sperm|germline|meiotic|meiosis|oocyt|oogen|vitell|egg antigen|CPEB",
  cycling_replication = "cell cycle|cyclin|mitotic|DNA replication|histone"
)
gene_sets <- lapply(program_patterns, function(pattern) {
  intersect(unique(functions$gene_id[grepl(pattern, functions$annotation_text,
                                            ignore.case = TRUE)]), rownames(obj))
})
gene_sets <- gene_sets[lengths(gene_sets) >= 2]

if (length(gene_sets) > 0) {
  set.seed(2026)
  module_obj <- AddModuleScore(obj, features = unname(gene_sets),
                               name = "validation_program_", seed = 2026)
  score_columns <- paste0("validation_program_", seq_along(gene_sets))
  names(score_columns) <- names(gene_sets)
  score_rows <- do.call(rbind, lapply(cluster_levels, function(cluster) {
    cells <- WhichCells(module_obj, idents = cluster)
    do.call(rbind, lapply(names(score_columns), function(program) {
      values <- module_obj[[score_columns[[program]]]][cells, 1]
      data.frame(cluster = cluster, cells = length(cells), module = program,
                 mean_score = mean(values), median_score = median(values))
    }))
  }))
  write.csv(score_rows, file.path(outdir, "all_cluster_functional_module_scores.csv"),
            row.names = FALSE)
}

writeLines(c(
  paste("Input RDS:", normalizePath(rds_file)),
  paste("Functional table:", normalizePath(master_file)),
  paste("Cells:", ncol(obj)),
  paste("Clusters:", length(unique(obj$seurat_clusters))),
  "", capture.output(sessionInfo())
), file.path(outdir, "manual_validation_provenance.txt"))

message("Manual validation completed")
message("Results: ", outdir)
