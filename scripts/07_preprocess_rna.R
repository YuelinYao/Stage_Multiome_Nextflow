#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
cfg <- read_pipeline_config(args[1], process = "RNA_processing")
settings <- cfg$RNA_processing
items <- read_stage_objects(args[-1], "RNA")
set.seed(123)
if (settings$integrate_rna && length(items) < 2L)
  stop_bad("RNA integration requires at least two samples; set integrate_rna=false for one sample")
# Preprocess each QC-passing sample in memory; only stage objects are saved.
items <- lapply(items, function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  if (min(ncol(obj), length(VariableFeatures(obj))) <= max(settings$RNA_dims))
    stop_bad("Insufficient cells/features for individual RNA PCA")
  obj <- RunPCA(obj, npcs = max(settings$RNA_dims), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = settings$RNA_dims, verbose = FALSE)
  obj <- FindClusters(obj, algorithm = settings$algorithm, resolution = settings$resolution,
    cluster.name = "sample_rna_clusters", verbose = FALSE)
  # Cell-cycle scores on the final, re-normalised cells; sex_prediction and
  # CellXY_db come from step 06 in the sample metadata.
  obj <- score_sample_cell_cycle(obj)
  SCTransform(obj, assay = "RNA", new.assay.name = "SCT",
    vst.flavor = "v2", method = if (requireNamespace("glmGamPoi", quietly = TRUE))
      "glmGamPoi" else "poisson", variable.features.n = 2000, verbose = FALSE)
})
annotations <- do.call(rbind, lapply(items, function(obj) {
  data.frame(cell_id = colnames(obj), obj[[]][, c("sample_id", "S.Score", "G2M.Score",
    "Phase", "sex_prediction", "CellXY_db"), drop = FALSE], row.names = NULL)
}))
write_tsv(annotations, "rna_cell_annotations.tsv")
features <- SelectIntegrationFeatures(items, nfeatures = 2000,
  assay = rep("SCT", length(items)))
items <- PrepSCTIntegration(items, anchor.features = features, verbose = FALSE)
obj <- if (length(items) == 1L) items[[1]] else
  merge(items[[1]], y = items[-1], merge.data = TRUE, merge.dr = FALSE)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
VariableFeatures(obj) <- features
required_pcs <- if (!is.null(cfg$component_count)) cfg$component_count else
  max(c(settings$RNA_dims, settings$RNA_umap_dims, cfg$WNN$RNA_dims))
if (min(ncol(obj), length(features)) <= required_pcs)
  stop_bad("Need more than ", required_pcs, " cells and SCT features for RNA PCA")
obj <- RunPCA(obj, features = features, npcs = required_pcs,
  reduction.name = "rna_merged_pca", reduction.key = "rnaPC_", verbose = FALSE)
obj <- cluster_qc_reduction(obj, "rna_merged_pca", settings$RNA_dims,
  settings$RNA_umap_dims, "rna_unintegrated", settings)
saveRDS(obj, "rna_unintegrated.rds")
if (settings$integrate_rna) {
  obj <- IntegrateLayers(obj, method = CCAIntegration, assay = "SCT",
    normalization.method = "SCT", features = features,
    orig.reduction = "rna_merged_pca", new.reduction = "integrated.cca",
    dims = settings$RNA_dims, dims.to.integrate = seq_len(required_pcs),
    k.weight = min(100L, min(vapply(items, ncol, integer(1))) - 1L), verbose = FALSE)
  obj <- cluster_qc_reduction(obj, "integrated.cca", settings$RNA_dims,
    settings$RNA_umap_dims, "rna_integrated", settings)
  obj <- PrepSCTFindMarkers(obj, verbose = FALSE)
  saveRDS(obj, "rna_integrated.rds")
}
