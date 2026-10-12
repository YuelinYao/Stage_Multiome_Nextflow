# Validate the 05 -> 06 handoff using real Seurat objects and saved QC tables.
# No Signac or DoubletFinder required. Run from the repository root:
# Rscript tests/reviewed_clusters.R
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error=identity), 'error'))
run_checks <- function() {
  file <- tempfile(fileext='.tsv'); on.exit(unlink(file), add=TRUE)
  cells <- sprintf('cell%03d', 1:120)
  tab <- data.frame(sample_id='WT_1', barcode=cells)
  for (metric in env$qc_metric_names) tab[[metric]] <- seq_along(cells)
  review <- tab; review$pre_doublet_cluster <- c(rep('0', 100), rep('7', 20))
  # Shuffle the reviewed table: assignment must use barcodes, never row order.
  env$write_tsv(review[120:1, ], file)
  labels <- env$read_reviewed_qc_clusters(file, 'WT_1', tab, cells)
  stopifnot(identical(names(labels), cells), identical(unname(labels), review$pre_doublet_cluster))
  obj <- CreateSeuratObject(matrix(1, 3, 120,
    dimnames=list(c('gene1','gene2','gene3'), cells)))
  obj$reviewed_qc_cluster <- unname(labels[colnames(obj)])
  # A different, unrelated clustering must never determine reviewed exclusions.
  obj$pre_doublet_clusters <- 'other_cluster'
  filtered <- env$remove_qc_clusters(obj, '7')
  stopifnot(ncol(filtered)==100, identical(colnames(filtered), cells[1:100]),
            all(filtered$reviewed_qc_cluster=='0'))
  expect_error(env$remove_qc_clusters(obj, 'unknown'))
  expect_error(env$read_reviewed_qc_clusters(file, 'WT_2', tab, cells))
  expect_error(env$read_reviewed_qc_clusters(file, 'WT_1', tab, cells[-1]))
  stale <- tab; stale$nCount_RNA[1] <- 999
  expect_error(env$read_reviewed_qc_clusters(file, 'WT_1', stale, cells))
  malformed <- review; malformed$pre_doublet_cluster[1] <- NA
  env$write_tsv(malformed, file)
  expect_error(env$read_reviewed_qc_clusters(file, 'WT_1', tab, cells))
  env$write_tsv(rbind(review, review[1, ]), file)
  expect_error(env$read_reviewed_qc_clusters(file, 'WT_1', tab, cells))
}
run_checks()
cat('All reviewed-cluster handoff checks passed.\n')
