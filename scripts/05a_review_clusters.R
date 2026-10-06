#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
sample <- args[1]; cfg <- read_analysis_config(args[4], sample)
assays <- make_assays(args[2], args[3])
tab <- qc_table(sample, assays)
keep <- apply_qc(tab, cfg$thresholds)
if (length(keep) < 100L) stop_bad("Fewer than 100 joint QC cells: ", sample)
rna <- preprocess_qc_rna(subset(assays$rna, cells = keep), cfg$clustering)
tab <- tab[match(colnames(rna), tab$barcode), , drop = FALSE]
tab$pre_doublet_cluster <- as.character(rna$pre_doublet_clusters)
write_tsv(tab, paste0(sample, ".cluster_cells.tsv"))
summary <- do.call(rbind, lapply(split(tab, tab$pre_doublet_cluster), function(x) {
  do.call(rbind, lapply(qc_metric_names, function(metric) {
    v <- x[[metric]]; v <- v[is.finite(v)]
    data.frame(sample_id = sample, cluster = x$pre_doublet_cluster[1], metric = metric,
      n = length(v), mean = mean(v), median = if (length(v)) median(v) else NA_real_)
  }))
}))
write_tsv(summary, paste0(sample, ".cluster_summary.tsv"))
plot_tab <- tab; plot_tab$sample_id <- paste0("Cluster ", tab$pre_doublet_cluster)
plot_qc_violins(plot_tab, paste0(sample, ".cluster_qc.pdf"), sample)
Idents(rna) <- "pre_doublet_clusters"
markers <- if (length(levels(Idents(rna))) > 1L)
  FindAllMarkers(rna, only.pos = TRUE, assay = "RNA", min.pct = 0.1,
    logfc.threshold = 0.25, verbose = FALSE) else data.frame()
# Apply significance to adjusted p-values explicitly.
if (nrow(markers)) {
  markers <- markers[is.finite(markers$p_val_adj) & markers$p_val_adj < 0.05,
                     , drop = FALSE]
}
write_tsv(markers, paste0(sample, ".cluster_markers.tsv"))
