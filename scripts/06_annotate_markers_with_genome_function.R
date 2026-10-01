#!/usr/bin/env Rscript

# Join Seurat marker genes to the final Swiss-Prot/eggNOG/InterPro master table.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 06_annotate_markers_with_genome_function.R MARKERS_CSV FUNCTION_MASTER_TSV OUTPUT_TSV")
}

marker_file <- args[1]
function_file <- args[2]
output_file <- args[3]
if (!file.exists(marker_file)) stop("Missing marker file: ", marker_file)
if (!file.exists(function_file)) stop("Missing function file: ", function_file)

markers <- read.csv(marker_file, check.names = FALSE, stringsAsFactors = FALSE)
functions <- read.delim(
  function_file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)
if (!"gene" %in% colnames(markers)) stop("Marker table lacks gene")
if (!"gene_id" %in% colnames(functions)) stop("Function table lacks gene_id")

if ("evidence_source_count" %in% colnames(functions)) {
  functions <- functions[order(functions$gene_id, -functions$evidence_source_count), ]
} else {
  functions <- functions[order(functions$gene_id), ]
}
functions <- functions[!duplicated(functions$gene_id), , drop = FALSE]

markers$.original_order <- seq_len(nrow(markers))
annotated <- merge(
  markers,
  functions,
  by.x = "gene",
  by.y = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
annotated <- annotated[order(annotated$.original_order), , drop = FALSE]
annotated$.original_order <- NULL

write.table(
  annotated,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

matched <- if ("protein_id" %in% colnames(annotated)) {
  sum(!is.na(annotated$protein_id) & nzchar(annotated$protein_id))
} else 0
cat("Marker rows:", nrow(markers), "\n")
cat("Rows with functional match:", matched, "\n")
cat("Percent matched:", round(100 * matched / nrow(markers), 2), "%\n")
cat("Output:", output_file, "\n")
