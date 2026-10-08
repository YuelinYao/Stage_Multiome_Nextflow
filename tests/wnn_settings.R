# Test the actual WNN stage with independent settings and synthetic Seurat
# reductions. No Signac, integration or real sequencing inputs required.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
# SCP is installed on the execution cluster, but absent from this test host.
# Substitute only its plotting entry points; Seurat WNN and PDF writing are real.
env$library <- function(package) stopifnot(as.character(substitute(package)) == 'SCP')
env$theme_blank <- function(xlab, ylab) ggplot2::theme_classic()
env$FeatureDimPlot <- function(srt, features, reduction, slot, pt.size, theme_use,
                              compare_features, color_blend_mode, label, label_insitu,
                              xlab, ylab, seed) {
  stopifnot(DefaultAssay(srt) == 'SCT', slot == 'data', reduction == 'wnn.umap',
            compare_features, color_blend_mode == 'blend', label, label_insitu,
            identical(features, c('Pou3f1', 'Foxa1', 'Pdgfra')))
  Seurat::FeaturePlot(srt, features = features[1], reduction = reduction, pt.size = pt.size)
}
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit(unlink(task_dir, recursive = TRUE), add = TRUE)
  set.seed(71)
  cells <- sprintf('cell%03d', 1:240)
  counts <- matrix(rpois(100 * 240, 3), 100, 240,
    dimnames = list(paste0('feature', 1:100), cells))
  rownames(counts)[1:3] <- c('Pou3f1', 'Foxa1', 'Pdgfra')
  rna <- CreateSeuratObject(counts, assay = 'RNA')
  rna[['SCT']] <- CreateAssayObject(data = log1p(counts))
  atac <- CreateSeuratObject(counts, assay = 'ATAC')
  rna$sample_id <- atac$sample_id <- rep(c('WT_1', 'WT_2'), each = 120)
  rna$barcode <- atac$barcode <- cells
  atac$nucleosome_signal <- 1; atac$TSS.enrichment <- 5
  reduction <- function(assay, key) {
    x <- matrix(rnorm(240 * 50), 240, 50,
      dimnames = list(cells, paste0(key, 1:50)))
    CreateDimReducObject(embeddings = x, key = key, assay = assay)
  }
  rna[['rna_merged_pca']] <- reduction('RNA', 'rnaPC_')
  rna[['integrated.cca']] <- reduction('RNA', 'cca_')
  atac[['lsi']] <- reduction('ATAC', 'LSI_')
  atac[['integrated_lsi']] <- reduction('ATAC', 'integratedLSI_')
  rna_file <- file.path(task_dir, 'rna.rds'); saveRDS(rna, rna_file)
  atac_file <- file.path(task_dir, 'atac.rds'); saveRDS(atac[, rev(cells)], atac_file)
  cfg <- list(RNA_processing = list(integrate_rna = TRUE, RNA_dims = 1:10),
              ATAC_processing = list(integrate_atac = FALSE, ATAC_dims = 2:20),
              WNN = list(RNA_dims = 1:35, ATAC_dims = 2:45,
                         algorithm = 1, resolution = 0.2, plot_point_size = 0.1,
                         lineage_assay = 'SCT',
                         lineage_markers = c('Pou3f1', 'Foxa1', 'Pdgfra')))
  env$read_pipeline_config <- function(...) cfg
  env$FindMultiModalNeighbors <- function(object, reduction.list, dims.list, ...) {
    expected <- list(if (cfg$RNA_processing$integrate_rna) 'integrated.cca' else 'rna_merged_pca',
                     if (cfg$ATAC_processing$integrate_atac) 'integrated_lsi' else 'lsi')
    stopifnot(identical(reduction.list, expected),
              identical(dims.list, list(1:35, 2:45)))
    Seurat::FindMultiModalNeighbors(object, reduction.list = reduction.list,
                                    dims.list = dims.list, ...)
  }
  env$FindClusters <- function(object, graph.name, resolution, algorithm, ...) {
    stopifnot(graph.name == 'wsnn', resolution == 0.2, algorithm == 1)
    Seurat::FindClusters(object, graph.name = graph.name, resolution = resolution,
                         algorithm = algorithm, ...)
  }
  env$args <- c(file.path(task_dir, 'wnn.rds'), 'json', rna_file, atac_file)
  for (expr in parse('scripts/12_wnn_stage.R')[-(1:2)]) eval(expr, env)
  result <- readRDS(env$args[1])
  stopifnot(DefaultAssay(result) == 'RNA',
    identical(as.character(result$orig.ident), as.character(rna$orig.ident)),
    all(file.info(file.path(task_dir, c('wnn_clusters_umap.pdf',
                                      'Cell_lineage_umap.pdf')))$size > 1000))
  stopifnot(identical(colnames(result), cells), 'wsnn_clusters' %in% colnames(result[[]]),
            identical(result@misc$wnn_reductions, c(RNA = 'integrated.cca', ATAC = 'lsi')))
  # Cover the other integration-switch combinations; WNN must select the
  # appropriate reduction from whichever saved modality objects it receives.
  for (switches in list(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE))) {
    cfg$RNA_processing$integrate_rna <- switches[1]
    cfg$ATAC_processing$integrate_atac <- switches[2]
    for (expr in parse('scripts/12_wnn_stage.R')[-(1:2)]) eval(expr, env)
    result <- readRDS(env$args[1])
    stopifnot(identical(result@misc$wnn_reductions,
      c(RNA = if (switches[1]) 'integrated.cca' else 'rna_merged_pca',
        ATAC = if (switches[2]) 'integrated_lsi' else 'lsi')))
  }
  # Missing gene symbols produce a report, preserving the analysis object.
  cfg$WNN$lineage_markers <- 'absent_gene'
  suppressWarnings(for (expr in parse('scripts/12_wnn_stage.R')[-(1:2)]) eval(expr, env))
  stopifnot(DefaultAssay(readRDS(env$args[1])) == 'RNA',
            file.info(file.path(task_dir, 'Cell_lineage_umap.pdf'))$size > 1000)
  cfg$WNN$lineage_markers <- c('Pou3f1', 'Foxa1', 'Pdgfra')
  # Different source cell order must be aligned; contradictory original
  # barcodes must fail before combining assays or computing neighbours.
  bad <- atac; bad$barcode[1] <- 'wrong_barcode'; saveRDS(bad, atac_file)
  expect_error <- function(expr) {
    stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
  }
  expect_error(for (expr in parse('scripts/12_wnn_stage.R')[-(1:2)]) eval(expr, env))
  saveRDS(atac[, rev(cells)], atac_file)
  # The standalone legacy combiner also retains the corrected ATAC reduction.
  env$args <- c(rna_file, atac_file, file.path(task_dir, 'combined.rds'))
  for (expr in parse('scripts/11_combine_modalities.R')[-(1:2)]) eval(expr, env)
  combined <- readRDS(env$args[3])
  stopifnot(identical(colnames(combined), cells),
            'integrated_lsi' %in% Reductions(combined))
}
run_checks()
cat('All independent WNN settings checks passed.\n')
