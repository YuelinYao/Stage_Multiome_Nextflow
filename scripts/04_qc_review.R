#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
# QC review writes measurements and plots; the user maintains the pipeline JSON.
stopifnot(length(args) >= 4L)
summary_file <- args[1]; candidates_file <- args[2]; plot_file <- args[3]
files <- args[-(1:3)]
message("04_qc_review.R arguments:\n",
        "  summary_file: ", summary_file, "\n",
        "  candidates_file: ", candidates_file, "\n",
        "  plot_file: ", plot_file, "\n",
        "  input QC TSV files:\n    ", paste(files, collapse = "\n    "))
tabs <- lapply(files, function(f) {
  # Preserve identifiers literally (including numeric-looking barcodes).
  tab <- read.delim(f, check.names = FALSE,
                   colClasses = c(sample_id = "character", barcode = "character"))
  required <- c("sample_id", "barcode", qc_metric_names)
  if (!all(required %in% names(tab)) || anyDuplicated(names(tab)))
    stop_bad("Missing or duplicated QC columns: ", f)
  if (!nrow(tab) || anyNA(tab[, c("sample_id", "barcode")]) ||
      any(!nzchar(tab$sample_id)) || any(!nzchar(tab$barcode)))
    stop_bad("Empty table or missing sample/barcode identifiers: ", f)
  if (length(unique(tab$sample_id)) != 1L || tab$sample_id[1] == "ALL_STAGE")
    stop_bad("Expected one sample per QC file; ALL_STAGE is reserved: ", f)
  if (!all(vapply(tab[qc_metric_names], is.numeric, logical(1))))
    stop_bad("QC metrics must be numeric: ", f)
  tab[, required, drop = FALSE]
})
# Nextflow tasks can finish in any order; keep combined tables/plots consistent.
tabs <- tabs[order(vapply(tabs, function(x) x$sample_id[1], character(1)))]
all <- do.call(rbind, tabs)
if (anyDuplicated(all[, c("sample_id", "barcode")])) stop_bad("Duplicate sample/barcode in QC files")
write_tsv(all, "all_cells_qc_metrics.tsv")
groups <- c(split(all, all$sample_id), list(ALL_STAGE = all))
summaries <- list(); candidates <- list()
for (sample in names(groups)) for (metric in qc_metric_names) {
  v <- groups[[sample]][[metric]]
  v <- v[is.finite(v)]
  if (!length(v)) stop_bad("No finite values: ", sample, "/", metric)
  med <- median(v); scale <- mad(v, constant = 1.4826)
  summaries[[length(summaries) + 1L]] <- data.frame(sample_id = sample, metric = metric,
    n = length(v), mean = mean(v), sd = sd(v), min = min(v), q01 = unname(quantile(v, .01)),
    q05 = unname(quantile(v, .05)), q25 = unname(quantile(v, .25)), median = med,
    q75 = unname(quantile(v, .75)), q95 = unname(quantile(v, .95)),
    q99 = unname(quantile(v, .99)), max = max(v), mad_scaled = scale,
    mad_3_lower = med - 3 * scale, mad_3_upper = med + 3 * scale,
    mad_5_lower = med - 5 * scale, mad_5_upper = med + 5 * scale,
    mad_7_lower = med - 7 * scale, mad_7_upper = med + 7 * scale)
  for (k in c(3L, 5L, 7L)) candidates[[length(candidates) + 1L]] <-
    data.frame(sample_id = sample, metric = metric, mad_multiplier = k,
               lower = med - k * scale, upper = med + k * scale)
}
write_tsv(do.call(rbind, summaries), summary_file)
write_tsv(do.call(rbind, candidates), candidates_file)
plot_qc_violins(all, plot_file, "QC distributions by sample")
