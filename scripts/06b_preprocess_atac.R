#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
cfg <- read_pipeline_config(args[1])
settings <- cfg$ATAC_processing
items <- read_stage_objects(args[-1], "ATAC")
if (!all(vapply(items, function(x) identical(rownames(x[["ATAC"]]),
    rownames(items[[1]][["ATAC"]])), logical(1))))
  stop_bad("ATAC inputs must share the same stage consensus peaks in the same order")
set.seed(123)
if (settings$integrate_atac && length(items) < 2L)
  stop_bad("ATAC integration requires at least two samples; set integrate_atac=false for one sample")
# Reciprocal-LSI anchors need sample-specific LSI, kept only in memory.
required_lsi <- max(c(settings$ATAC_dims, settings$ATAC_integrate_dims, cfg$WNN$ATAC_dims))
items <- lapply(items, function(obj) {
  DefaultAssay(obj) <- "ATAC"
  obj <- FindTopFeatures(obj, min.cutoff = settings$min_cutoff)
  obj <- RunTFIDF(obj, verbose = FALSE)
  if (min(ncol(obj), length(VariableFeatures(obj))) <= required_lsi)
    stop_bad("Need more than ", required_lsi, " cells and variable peaks for sample LSI")
  RunSVD(obj, n = required_lsi, verbose = FALSE)
})
obj <- if (length(items) == 1L) items[[1]] else
  merge(items[[1]], y = items[-1], merge.data = TRUE, merge.dr = FALSE)
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = settings$min_cutoff)
obj <- RunTFIDF(obj, verbose = FALSE)
required_lsi <- max(c(settings$ATAC_dims, settings$ATAC_integrate_dims, cfg$WNN$ATAC_dims))
if (min(ncol(obj), length(VariableFeatures(obj))) <= required_lsi)
  stop_bad("Need more than ", required_lsi, " cells and variable peaks for LSI")
obj <- RunSVD(obj, n = required_lsi, verbose = FALSE)
obj <- cluster_qc_reduction(obj, "lsi", settings$ATAC_dims,
  settings$ATAC_dims, "atac_unintegrated", settings)
saveRDS(obj, "atac_unintegrated.rds")
if (settings$integrate_atac) {
  anchors <- FindIntegrationAnchors(object.list = items,
    reduction = "rlsi", dims = settings$ATAC_dims)
  # IntegrateEmbeddings returns a new merged object; copy its corrected reduction
  # to retain the unintegrated graphs, clusters, UMAP and original merged assays.
  corrected <- IntegrateEmbeddings(anchorset = anchors, reductions = obj[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = settings$ATAC_integrate_dims,
    k.weight = min(100L, min(vapply(items, ncol, integer(1))) - 1L))
  obj[["integrated_lsi"]] <- corrected[["integrated_lsi"]]
  obj <- cluster_qc_reduction(obj, "integrated_lsi", settings$ATAC_dims,
    settings$ATAC_dims, "atac_integrated", settings)
  saveRDS(obj, "atac_integrated.rds")
}
saveRDS(obj, "atac_stage.rds")
