# Exercise the actual 07 script, including SCT and CCA integration, with
# synthetic RNA counts. Only cell-cycle/sex annotation is stubbed because
# cellXY is unavailable locally. Run from the repository root.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  original_dir <- getwd()
  script <- normalizePath('scripts/07_preprocess_rna.R')
  cfg <- jsonlite::fromJSON('examples/pipeline_settings.example.json', simplifyVector = FALSE)
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit({ setwd(original_dir); unlink(task_dir, recursive = TRUE) }, add = TRUE)
  set.seed(29)
  # Shared biological groups with modest sample-specific expression effects.
  gene_mean <- rgamma(900, shape = 2, rate = 0.5)
  files <- vapply(c('A', 'B'), function(id) {
    cells <- paste0(id, '_', seq_len(140))
    means <- outer(gene_mean, runif(140, 0.8, 1.2))
    means[1:100, 1:70] <- means[1:100, 1:70] * 3
    if (id == 'B') means[101:150, ] <- means[101:150, ] * 1.2
    counts <- matrix(rpois(length(means), means), 900, 140,
      dimnames = list(paste0('gene', seq_len(900)), cells))
    obj <- CreateSeuratObject(counts, assay = 'RNA'); obj$sample_id <- id
    # Step 06 adds the cellXY columns before step 07 reads the object.
    obj$sex_prediction <- 'Male'; obj$CellXY_db <- 'Singlet'
    file <- file.path(task_dir, paste0(id, '.rds')); saveRDS(obj, file); file
  }, character(1))
  cfg$RNA_processing$algorithm <- 1
  cfg$RNA_processing$resolution <- 0.2
  cfg$RNA_processing$RNA_dims <- '1:10'
  cfg$RNA_processing$RNA_umap_dims <- '1:20'
  cfg$WNN$RNA_dims <- '1:30'
  cfg <- list(process = 'RNA_processing', RNA_processing = cfg$RNA_processing, component_count = 30)
  config <- file.path(task_dir, 'settings.json')
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, null = 'null')
  # Step 07 scores the cell cycle on each sample; cellXY columns come from step 06.
  annotated <- 0L
  env$score_sample_cell_cycle <- function(obj) {
    annotated <<- annotated + 1L
    obj$S.Score <- 1; obj$G2M.Score <- 1; obj$Phase <- 'S'; obj
  }
  env$annotate_cellxy <- function(obj) stop('cellXY runs in step 06 only')
  env$args <- c(config, files)
  setwd(task_dir)
  for (expr in parse(script)[-(1:2)]) eval(expr, env)
  result <- readRDS('rna_integrated.rds')
  stopifnot(ncol(result) == 280,
            ncol(Embeddings(result, 'rna_merged_pca')) == 30,
            ncol(Embeddings(result, 'integrated.cca')) == 30,
            all(c('rna_unintegrated_umap', 'rna_integrated_umap') %in% Reductions(result)),
            all(c('rna_unintegrated_clusters', 'rna_integrated_clusters') %in% colnames(result[[]])),
            file.exists('rna_unintegrated.rds'), file.exists('rna_integrated.rds'),
            !file.exists('rna_stage.rds'))
  annotations <- read.delim('rna_cell_annotations.tsv')
  stopifnot(nrow(annotations) == 280, !anyDuplicated(annotations$cell_id),
            all(c('sex_prediction', 'CellXY_db') %in% names(annotations)),
            all(annotations$CellXY_db == 'Singlet'), annotated == 2L,
            all(annotations$Phase == 'S'), all(annotations$sex_prediction == 'Male'))
  # A fresh single-sample run without integration saves only the unintegrated
  # object; it is also the object that Nextflow supplies to WNN.
  unlink(c('rna_unintegrated.rds', 'rna_integrated.rds', 'rna_cell_annotations.tsv'))
  cfg$RNA_processing$integrate_rna <- FALSE
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, null = 'null')
  env$args <- c(config, files[1])
  for (expr in parse(script)[-(1:2)]) eval(expr, env)
  result <- readRDS('rna_unintegrated.rds')
  stopifnot(ncol(result) == 140, 'rna_merged_pca' %in% Reductions(result),
            !('integrated.cca' %in% Reductions(result)),
            !file.exists('rna_integrated.rds'), !file.exists('rna_stage.rds'))
}
run_checks()
cat('All RNA SCT/CCA integration checks passed.\n')
