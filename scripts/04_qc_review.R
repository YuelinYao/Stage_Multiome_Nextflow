#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5L)
summary_file <- args[1]; candidates_file <- args[2]; json_file <- args[3]
plot_file <- args[4]
files <- args[-(1:4)]
tabs <- lapply(files, function(f) read.delim(f, check.names = FALSE))
all <- do.call(rbind, tabs)
if (anyDuplicated(paste(all$sample_id, all$barcode))) stop_bad("Duplicate sample/barcode in QC files")
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
template <- list(
  note = "Fill the min/max values after reviewing plots and MAD candidates. null disables one side. Do not use this template unchanged.",
  thresholds = setNames(lapply(qc_metric_names, function(nm) list(min = NULL, max = NULL)), qc_metric_names),
  sample_overrides = list(),
  doublet = list(expected_rate = 0.05, sample_rates = list(), exclude_clusters = list()),
  clustering = list(algorithm = 3L, resolution = 2,
                    RNA_dims = "1:20", RNA_umap_dims = "1:50",
                    ATAC_dims = "2:40", ATAC_integrate_dims = "1:40",
                    integrate_rna = TRUE, integrate_atac = TRUE, min_cutoff = 10L))
jsonlite::write_json(template, json_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
