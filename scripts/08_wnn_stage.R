#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
outfile <- args[1]
cfg <- read_pipeline_config(args[2])
settings <- cfg$WNN
obj <- readRDS(args[3]); atac <- readRDS(args[4])
if (!setequal(colnames(obj), colnames(atac))) stop_bad("Stage RNA/ATAC cells differ")
atac <- atac[, colnames(obj)]
if (!identical(as.character(obj$sample_id), as.character(atac$sample_id)))
  stop_bad("Stage RNA/ATAC sample assignments differ")
obj[["ATAC"]] <- atac[["ATAC"]]
for (key in c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"))
  obj[[key]] <- atac[[key]][colnames(obj), 1, drop = TRUE]
for (reduction in Reductions(atac)) obj[[reduction]] <- atac[[reduction]]
for (key in c("atac_unintegrated_clusters", "atac_integrated_clusters"))
  if (key %in% colnames(atac[[]])) obj[[key]] <- atac[[key]][, 1]
rna_reduction <- if (cfg$RNA_processing$integrate_rna) "integrated.cca" else "rna_merged_pca"
atac_reduction <- if (cfg$ATAC_processing$integrate_atac) "integrated_lsi" else "lsi"
for (reduction in c(rna_reduction, atac_reduction))
  if (!(reduction %in% Reductions(obj))) stop_bad("Missing stage reduction: ", reduction)
if (ncol(Embeddings(obj, rna_reduction)) < max(settings$RNA_dims) ||
    ncol(Embeddings(obj, atac_reduction)) < max(settings$ATAC_dims))
  stop_bad("Insufficient dimensions for WNN")
set.seed(123)
DefaultAssay(obj) <- "RNA"
obj <- FindMultiModalNeighbors(obj,
  reduction.list = list(rna_reduction, atac_reduction),
  dims.list = list(settings$RNA_dims, settings$ATAC_dims),
  modality.weight.name = c("RNA.weight", "ATAC.weight"))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", resolution = settings$resolution,
  cluster.name = "wsnn_clusters", algorithm = settings$algorithm)
obj@misc$wnn_reductions <- c(RNA = rna_reduction, ATAC = atac_reduction)
saveRDS(obj, outfile)
