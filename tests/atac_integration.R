# Test actual 08 routing and real Seurat reciprocal-LSI/embedding integration.
# Signac preprocessing is replaced with synthetic feature/LSI preparation;
# this does not validate ChromatinAssay fragment handling or Signac itself.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  original_dir <- getwd(); script <- normalizePath('scripts/08_preprocess_atac.R')
  cfg <- jsonlite::fromJSON('examples/pipeline_settings.example.json', simplifyVector = FALSE)
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit({ setwd(original_dir); unlink(task_dir, recursive = TRUE) }, add = TRUE)
  previous <- options(Seurat.object.assay.version = 'v3')
  on.exit(options(previous), add = TRUE)
  set.seed(91)
  features <- paste0('chr1-', seq(100, by = 100, length.out = 320), '-',
                     seq(150, by = 100, length.out = 320))
  files <- vapply(c('A', 'B'), function(id) {
    counts <- matrix(rpois(320 * 60, 2), 320, 60,
      dimnames = list(features, paste0(id, '_', seq_len(60))))
    counts[1:50, 1:30] <- counts[1:50, 1:30] + 3
    obj <- CreateSeuratObject(counts, assay = 'ATAC'); obj$sample_id <- id
    file <- file.path(task_dir, paste0(id, '.rds')); saveRDS(obj, file); file
  }, character(1))
  cfg$ATAC_processing$ATAC_dims <- '2:10'
  cfg$ATAC_processing$ATAC_integrate_dims <- '1:80'
  cfg$ATAC_processing$algorithm <- 1
  cfg$ATAC_processing$resolution <- 0.2
  cfg$WNN$ATAC_dims <- '2:80'
  cfg <- list(process = 'ATAC_processing', ATAC_processing = cfg$ATAC_processing, component_count = 80)
  config <- file.path(task_dir, 'settings.json')
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, null = 'null')
  env$FindTopFeatures <- function(object, min.cutoff) {
    stopifnot(identical(min.cutoff, 'q5'))
    # Fixture features are supplied to isolate script routing and integration.
    VariableFeatures(object) <- rownames(object); object
  }
  env$RunTFIDF <- function(object, verbose) {
    # Synthetic numeric data for reciprocal projections, not Signac TF-IDF.
    NormalizeData(object, verbose = FALSE)
  }
  requested <- integer()
  env$RunSVD <- function(object, n, verbose) {
    requested <<- c(requested, n)
    data <- as.matrix(GetAssayData(object, assay = 'ATAC', layer = 'data'))
    decomposition <- svd(t(data), nu = n, nv = n)
    embeddings <- sweep(decomposition$u[, seq_len(n), drop = FALSE], 2,
                         decomposition$d[seq_len(n)], '*')
    rownames(embeddings) <- colnames(object); colnames(embeddings) <- paste0('LSI_', seq_len(n))
    loadings <- decomposition$v[, seq_len(n), drop = FALSE]
    rownames(loadings) <- rownames(object); colnames(loadings) <- colnames(embeddings)
    object[['lsi']] <- CreateDimReducObject(embeddings = embeddings, loadings = loadings,
      key = 'LSI_', assay = 'ATAC',
      misc = list(d = decomposition$d[seq_len(n)],
                  u = decomposition$u[, seq_len(n), drop = FALSE]))
    object
  }
  # Motif/gene-activity assays need real peaks, fragments and BSgenome; add stand-in
  # assays here and check that both saved objects keep them.
  feature_calls <- 0L
  env$add_atac_feature_assays <- function(obj, cpus) {
    stopifnot(cpus >= 1L); feature_calls <<- feature_calls + 1L
    fake <- matrix(1, 3, ncol(obj), dimnames = list(c('MA1', 'MA2', 'MA3'), colnames(obj)))
    obj[['chromvar']] <- CreateAssayObject(data = fake)
    obj[['chromvar_rename']] <- CreateAssayObject(data = fake)
    obj[['gene.activities']] <- CreateAssayObject(counts = fake)
    obj
  }
  feature_assays <- c('chromvar', 'chromvar_rename', 'gene.activities')
  env$args <- c(config, files)
  setwd(task_dir)
  expressions <- parse(script)[-(1:2)]
  for (expr in expressions) eval(expr, env)
  result <- readRDS('atac_integrated.rds')
  stopifnot(identical(requested, c(10L, 10L, 80L)), ncol(result) == 120,
            ncol(Embeddings(result, 'lsi')) == 80,
            ncol(Embeddings(result, 'integrated_lsi')) == 80,
            identical(rownames(Embeddings(result, 'integrated_lsi')), colnames(result)),
            result@misc$atac_integration_k_weight < 59,
            identical(rownames(result[['ATAC']]), features),
            all(c('atac_unintegrated_umap', 'atac_integrated_umap') %in% Reductions(result)),
            file.exists('atac_integrated.rds'), !file.exists('atac_stage.rds'),
            feature_calls == 1L, all(feature_assays %in% Assays(result)),
            all(feature_assays %in% Assays(readRDS('atac_unintegrated.rds'))),
            DefaultAssay(result) == 'ATAC')
  original_counts <- do.call(cbind, lapply(files, function(file)
    GetAssayData(readRDS(file), assay = 'ATAC', layer = 'counts')))
  stopifnot(identical(GetAssayData(result, assay = 'ATAC', layer = 'counts'), original_counts))
  # Without integration, sample LSI and the unused integration range are skipped.
  unlink(c('atac_integrated.rds', 'atac_unintegrated.rds'))
  cfg$ATAC_processing$integrate_atac <- FALSE
  cfg$ATAC_processing$ATAC_integrate_dims <- NULL
  cfg$component_count <- 30
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, null = 'null')
  requested <- integer()
  for (expr in expressions) eval(expr, env)
  result <- readRDS('atac_unintegrated.rds')
  stopifnot(identical(requested, 30L), ncol(Embeddings(result, 'lsi')) == 30,
            !('integrated_lsi' %in% Reductions(result)),
            !file.exists('atac_integrated.rds'), !file.exists('atac_stage.rds'))
}
run_checks()
# chromVAR rows are renamed from motif ID to TF name; repeated names become unique.
deviations <- matrix(1:6, 3, 2, dimnames = list(c('MA1', 'MA2', 'MA3'), c('c1', 'c2')))
renamed <- env$rename_motif_rows(deviations, list(MA1 = 'Ctcf', MA2 = 'Foxa1', MA3 = 'Ctcf'))
stopifnot(identical(rownames(renamed), c('Ctcf', 'Foxa1', 'Ctcf.1')),
          identical(unname(renamed), unname(deviations)),
          inherits(tryCatch(env$rename_motif_rows(deviations, list(MA1 = 'Ctcf')),
                            error = identity), 'error'))
cat('All ATAC stage/integration checks passed.\n')
