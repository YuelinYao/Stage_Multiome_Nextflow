# Exercise the actual 06 orchestration with real Seurat preprocessing and
# mocked DoubletFinder calls; no Signac or DoubletFinder installation required.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit(unlink(task_dir, recursive = TRUE), add = TRUE)
  set.seed(19)
  cells <- sprintf('cell%03d', 1:140)
  counts <- matrix(rpois(300 * 140, 4), 300, 140,
    dimnames = list(paste0('gene', 1:300), cells))
  obj <- CreateSeuratObject(counts)
  tab <- data.frame(sample_id = 'WT_1', barcode = cells)
  for (metric in env$qc_metric_names) tab[[metric]] <- rep(1, 140)
  tab$nCount_RNA[1:10] <- 0
  limits <- setNames(lapply(env$qc_metric_names, function(x) list(min = 1, max = 2)),
                    env$qc_metric_names)
  review <- tab[11:140, ]; review$pre_doublet_cluster <- c(rep('7', 20), rep('0', 110))
  review_file <- file.path(task_dir, 'review.tsv'); env$write_tsv(review, review_file)
  env$read_pipeline_config <- function(...) list(thresholds = limits,
    review_QC_cluster = list(RNA_dims = 1:2, resolution = 9, algorithm = 3),
    doublet = list(PC_dims = 1:15, algorithm = 1, resolution = 0.1, expected_rate = 'auto', exclude_clusters = '7'))
  arc <- file.path(task_dir, 'outs'); dir.create(arc)
  write.csv(data.frame('Estimated number of cells' = '5,000', check.names = FALSE),
            file.path(arc, 'summary.csv'), row.names = FALSE)
  env$make_assays <- function(...) list(rna = obj, atac = obj)
  env$qc_table <- function(...) tab
  preprocess <- env$preprocess_qc_rna
  env$preprocess_qc_rna <- function(obj, settings, default_pca) {
    stopifnot(settings$algorithm == 1, settings$resolution == 0.1,
              identical(settings$RNA_dims, 1:15), default_pca,
              !('pre_doublet_clusters' %in% colnames(obj[[]])),
              all(obj$reviewed_qc_cluster == '0'))
    preprocess(obj, settings, default_pca)
  }
  env$paramSweep <- function(seu, PCs, sct) {
    stopifnot(identical(PCs, 1:15), identical(sct, FALSE), ncol(seu) == 110,
              ncol(Embeddings(seu, 'pca')) == 50,
              all(seu$reviewed_qc_cluster == '0'))
    list()
  }
  env$summarizeSweep <- function(...) data.frame()
  # Verify factor pK conversion and exclusion of a non-finite sweep result.
  env$find.pK <- function(...) data.frame(pK = factor(c('0.01', '0.02')),
                                        BCmetric = c(2, NA_real_))
  env$modelHomotypic <- function(annotations) {
    stopifnot(length(annotations) == 110); 0.2
  }
  env$doubletFinder <- function(seu, PCs, pK, nExp, sct) {
    stopifnot(identical(PCs, 1:15), pK == 0.01, nExp == 3, !sct)
    seu$DF.classifications_test <- c(rep('Singlet', 107), rep('Doublet', 3)); seu
  }
  env$args <- c('WT_1', arc, 'bed', 'json', file.path(task_dir, 'rna.rds'),
                file.path(task_dir, 'atac.rds'), file.path(task_dir, 'qc.tsv'), review_file)
  exprs <- parse('scripts/06_qc_doubletfinder.R')
  for (expr in exprs[-(1:2)]) {
    if (is.call(expr) && identical(expr[[1]], as.name('suppressPackageStartupMessages'))) next
    eval(expr, env)
  }
  report <- read.delim(env$args[7])
  stopifnot(report$raw_joint == 140, report$qc_pass == 130, report$qc_removed == 10,
            report$cluster_pass == 110, report$cluster_removed == 20,
            report$df_singlets == 107, report$df_removed == 3,
            report$arc_detected_cells == 5000, report$expected_rate == 0.039,
            report$expected_doublets_unadjusted == 4, report$expected_doublets == 3,
            report$df_pcs == paste(1:15, collapse = ','), report$pca_computed == 50)
  rna <- readRDS(env$args[5]); atac <- readRDS(env$args[6])
  barcodes <- read.delim(file.path(task_dir, 'WT_1.singlet_barcodes.tsv'),
                        colClasses = 'character')
  stopifnot(identical(names(barcodes), c('sample_id', 'barcode')),
            all(barcodes$sample_id == 'WT_1'), identical(barcodes$barcode, unname(rna$barcode)))
  stopifnot(file.info(file.path(task_dir, 'WT_1_Doublet.pdf'))$size > 0,
            all(rna$DoubletFinder == 'Singlet'))
  stopifnot(identical(colnames(rna), colnames(atac)), ncol(rna) == 107,
            all(rna$sample_id == 'WT_1'), identical(unname(rna$barcode), cells[31:137]))
}
run_checks()
cat('All DoubletFinder stage checks passed.\n')
