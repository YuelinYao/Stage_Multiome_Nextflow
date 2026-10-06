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
if (settings$integrate_atac) {
  sample_lsi <- max(settings$ATAC_dims)
  items <- lapply(items, function(obj) {
    DefaultAssay(obj) <- "ATAC"
    obj <- FindTopFeatures(obj, min.cutoff = settings$min_cutoff)
    obj <- RunTFIDF(obj, verbose = FALSE)
    if (min(ncol(obj), length(VariableFeatures(obj))) <= sample_lsi)
      stop_bad("Need more than ", sample_lsi, " cells and variable peaks for sample LSI")
    RunSVD(obj, n = sample_lsi, verbose = FALSE)
  })
}
obj <- if (length(items) == 1L) items[[1]] else
  merge(items[[1]], y = items[-1], merge.data = TRUE, merge.dr = FALSE)
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = settings$min_cutoff)
obj <- RunTFIDF(obj, verbose = FALSE)
required_lsi <- max(c(settings$ATAC_dims, cfg$WNN$ATAC_dims,
  if (settings$integrate_atac) settings$ATAC_integrate_dims))
if (min(ncol(obj), length(VariableFeatures(obj))) <= required_lsi)
  stop_bad("Need more than ", required_lsi, " cells and variable peaks for LSI")
obj <- RunSVD(obj, n = required_lsi, verbose = FALSE)
obj <- cluster_qc_reduction(obj, "lsi", settings$ATAC_dims,
  settings$ATAC_dims, "atac_unintegrated", settings)
saveRDS(obj, "atac_unintegrated.rds")
if (settings$integrate_atac) {
  # Signac tutorial: anchors use 2:30, correction retains 1:30, and downstream
  # UMAP uses 2:30. Retaining component 1 does not include it in neighbours or
  # UMAP when ATAC_dims starts at 2; it preserves component numbering.
  # https://stuartlab.org/signac/articles/integrate_atac
  anchors <- FindIntegrationAnchors(object.list = items,
    anchor.features = rownames(obj[["ATAC"]]),
    reduction = "rlsi", dims = settings$ATAC_dims, verbose = FALSE)
  # Weighting requires distinct anchor cells, not just a sufficient sample
  # size. Use a conservative bound valid for either direction of each pair.
  anchor_table <- as.data.frame(methods::slot(anchors, "anchors"))
  if (!nrow(anchor_table)) stop_bad("No reciprocal-LSI integration anchors found")
  pair <- interaction(anchor_table$dataset1, anchor_table$dataset2, drop = TRUE)
  anchor_cells <- vapply(split(anchor_table$cell2, pair),
    function(cells) length(unique(cells)), integer(1))
  if (min(anchor_cells) < 3L)
    stop_bad("Fewer than three distinct anchor cells in an ATAC sample pair")
  k_weight <- min(100L, min(vapply(items, ncol, integer(1))) - 1L,
                  min(anchor_cells) - 1L)
  message("ATAC integration uses k.weight = ", k_weight)
  # IntegrateEmbeddings returns a new merged object; copy its corrected reduction
  # to retain the unintegrated graphs, clusters, UMAP and original merged assays.
  corrected <- IntegrateEmbeddings(anchorset = anchors, reductions = obj[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = settings$ATAC_integrate_dims,
    k.weight = k_weight, verbose = FALSE)
  if (!setequal(rownames(Embeddings(corrected, "integrated_lsi")), colnames(obj)))
    stop_bad("Integrated ATAC embeddings do not match stage cells")
  obj[["integrated_lsi"]] <- corrected[["integrated_lsi"]]
  obj@misc$atac_integration_k_weight <- k_weight
  obj <- cluster_qc_reduction(obj, "integrated_lsi", settings$ATAC_dims,
    settings$ATAC_dims, "atac_integrated", settings)
  saveRDS(obj, "atac_integrated.rds")
}
