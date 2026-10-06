# Exercise the actual 06a script, including SCT and CCA integration, with
# synthetic RNA counts. Only cell-cycle/sex annotation is stubbed because
# cellXY is unavailable locally. Run from the repository root.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  original_dir <- getwd()
  script <- normalizePath('scripts/06a_preprocess_rna.R')
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
    file <- file.path(task_dir, paste0(id, '.rds')); saveRDS(obj, file); file
  }, character(1))
  cfg$RNA_processing$algorithm <- 1
  cfg$RNA_processing$resolution <- 0.2
  cfg$RNA_processing$RNA_dims <- '1:10'
  cfg$RNA_processing$RNA_umap_dims <- '1:20'
  cfg$WNN$RNA_dims <- '1:30'
  config <- file.path(task_dir, 'settings.json')
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, null = 'null')
  env$annotate_sample_rna <- function(obj) {
    obj$S.Score <- 0; obj$G2M.Score <- 0
    obj$Phase <- 'G1'; obj$sex_prediction <- NA_character_; obj
  }
  env$args <- c(config, files)
  setwd(task_dir)
  for (expr in parse(script)[-(1:2)]) eval(expr, env)
  result <- readRDS('rna_stage.rds')
  stopifnot(ncol(result) == 280,
            ncol(Embeddings(result, 'rna_merged_pca')) == 30,
            ncol(Embeddings(result, 'integrated.cca')) == 30,
            all(c('rna_unintegrated_umap', 'rna_integrated_umap') %in% Reductions(result)),
            all(c('rna_unintegrated_clusters', 'rna_integrated_clusters') %in% colnames(result[[]])),
            file.exists('rna_unintegrated.rds'), file.exists('rna_integrated.rds'))
  annotations <- read.delim('rna_cell_annotations.tsv')
  stopifnot(nrow(annotations) == 280, !anyDuplicated(annotations$cell_id))
}
run_checks()
cat('All RNA SCT/CCA integration checks passed.\n')
