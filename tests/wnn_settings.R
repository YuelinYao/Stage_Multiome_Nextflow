# Test the actual WNN stage with independent settings and synthetic Seurat
# reductions. No Signac, integration or real sequencing inputs required.
suppressPackageStartupMessages(library(Seurat))
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit(unlink(task_dir, recursive = TRUE), add = TRUE)
  set.seed(71)
  cells <- sprintf('cell%03d', 1:240)
  counts <- matrix(rpois(100 * 240, 3), 100, 240,
    dimnames = list(paste0('feature', 1:100), cells))
  rna <- CreateSeuratObject(counts, assay = 'RNA')
  atac <- CreateSeuratObject(counts, assay = 'ATAC')
  rna$sample_id <- atac$sample_id <- 'WT_1'
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
  atac_file <- file.path(task_dir, 'atac.rds'); saveRDS(atac, atac_file)
  cfg <- list(RNA_processing = list(integrate_rna = TRUE, RNA_dims = 1:10),
              ATAC_processing = list(integrate_atac = FALSE, ATAC_dims = 2:20),
              WNN = list(RNA_dims = 1:35, ATAC_dims = 2:45,
                         algorithm = 1, resolution = 0.2))
  env$read_pipeline_config <- function(...) cfg
  env$FindMultiModalNeighbors <- function(object, reduction.list, dims.list, ...) {
    stopifnot(identical(reduction.list, list('integrated.cca', 'lsi')),
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
  for (expr in parse('scripts/08_wnn_stage.R')[-(1:2)]) eval(expr, env)
  result <- readRDS(env$args[1])
  stopifnot(identical(colnames(result), cells), 'wsnn_clusters' %in% colnames(result[[]]),
            identical(result@misc$wnn_reductions, c(RNA = 'integrated.cca', ATAC = 'lsi')))
}
run_checks()
cat('All independent WNN settings checks passed.\n')
