# Exercise legacy modality combination and actual final RNA marker testing.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  original_dir <- getwd(); root <- tempfile(); dir.create(root)
  legacy <- normalizePath('scripts/11_combine_modalities.R')
  markers <- normalizePath('scripts/14_cluster_markers.R')
  on.exit({ setwd(original_dir); unlink(root, recursive = TRUE) }, add = TRUE)
  set.seed(22)
  counts <- matrix(rpois(80 * 120, 2), 80, 120,
    dimnames = list(paste0('gene', 1:80), paste0('cell', 1:120)))
  counts[1:5, 1:60] <- counts[1:5, 1:60] + 20
  counts[6:10, 61:120] <- counts[6:10, 61:120] + 20
  rna <- NormalizeData(CreateSeuratObject(counts), verbose = FALSE)
  atac <- CreateSeuratObject(counts, assay = 'ATAC')
  rna$sample_id <- atac$sample_id <- rep(c('WT_1', 'WT_2'), each = 60)
  rna$barcode <- atac$barcode <- paste0('barcode', 1:120)
  rna$wsnn_clusters <- rep(c('0', '1'), each = 60)
  atac$nucleosome_signal <- 1; atac$TSS.enrichment <- 5
  embeddings <- matrix(rnorm(120 * 10), 120, 10,
    dimnames = list(colnames(atac), paste0('LSI_', 1:10)))
  atac[['lsi']] <- CreateDimReducObject(embeddings, key = 'LSI_', assay = 'ATAC')
  atac[['integrated_lsi']] <- CreateDimReducObject(embeddings, key = 'iLSI_', assay = 'ATAC')
  rna_file <- file.path(root, 'rna.rds'); saveRDS(rna, rna_file)
  atac_file <- file.path(root, 'atac.rds'); saveRDS(atac[, rev(colnames(atac))], atac_file)
  combined <- file.path(root, 'combined.rds')
  env$args <- c(rna_file, atac_file, combined)
  expressions <- parse(legacy)[-(1:2)]
  for (expr in expressions) eval(expr, env)
  obj <- readRDS(combined)
  stopifnot(identical(colnames(obj), colnames(rna)), 'ATAC' %in% Assays(obj),
            all(c('lsi', 'integrated_lsi') %in% Reductions(obj)),
            identical(Embeddings(obj, 'lsi'), Embeddings(atac, 'lsi')))
  # Mismatched sample identity fails even when cell names match.
  bad <- atac; bad$sample_id <- 'wrong'; saveRDS(bad, atac_file)
  stopifnot(inherits(tryCatch({ for (expr in expressions) eval(expr, env) }, error = identity), 'error'))
  setwd(root)
  env$args <- combined
  for (expr in parse(markers)[-(1:2)]) eval(expr, env)
  result <- read.delim('cluster_markers.tsv')
  stopifnot(nrow(result) > 0, all(is.finite(result$p_val_adj) & result$p_val_adj < 0.05),
            all(paste0('gene', 1:10) %in% result$gene))
  obj$wsnn_clusters <- '0'; saveRDS(obj, combined)
  for (expr in parse(markers)[-(1:2)]) eval(expr, env)
  empty <- read.delim('cluster_markers.tsv')
  stopifnot(nrow(empty) == 0, all(c('p_val_adj', 'cluster', 'gene') %in% names(empty)))
}
run_checks()
cat('All legacy-modality and final-marker checks passed.\n')
