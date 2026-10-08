#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(SCP))
stopifnot(length(args) == 4L)
outfile <- args[1]
cfg <- read_pipeline_config(args[2], process = "WNN")
settings <- cfg$WNN
obj <- readRDS(args[3]); atac <- readRDS(args[4])
if (!setequal(colnames(obj), colnames(atac))) stop_bad("Stage RNA/ATAC cells differ")
atac <- atac[, colnames(obj)]
if (!identical(as.character(obj$sample_id), as.character(atac$sample_id)))
  stop_bad("Stage RNA/ATAC sample assignments differ")
if (is.null(obj$barcode) || is.null(atac$barcode) ||
    !identical(as.character(obj$barcode), as.character(atac$barcode)))
  stop_bad("Stage RNA/ATAC original barcodes differ or are missing")
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
# Keep the existing reduction name in saved objects and all plots.
reduction <- "wnn.umap"
plot_dir <- dirname(outfile)
sample_groups <- levels(factor(obj$sample_id))
cluster_groups <- levels(factor(obj$wsnn_clusters))
p_umap1 <- DimPlot(obj, reduction = reduction, group.by = "sample_id",
  pt.size = settings$plot_point_size, cols = plot_group_colors(sample_groups)) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1,
    override.aes = list(size = 6))) + ggplot2::ggtitle("") +
  theme_blank(xlab = "wnnUMAP_1", ylab = "wnnUMAP_2")
p_umap2 <- DimPlot(obj, reduction = reduction, group.by = "wsnn_clusters",
  label = TRUE, repel = TRUE, pt.size = settings$plot_point_size,
  cols = plot_group_colors(cluster_groups)) + NoLegend() + ggplot2::ggtitle("") +
  theme_blank(xlab = "wnnUMAP_1", ylab = "wnnUMAP_2")
ggplot2::ggsave(file.path(plot_dir, "wnn_clusters_umap.pdf"),
  plot = patchwork::wrap_plots(p_umap1, p_umap2, nrow = 1), width = 9, height = 4.2)

# SCT data contains normalized expression; scale.data contains residuals and
# may omit these genes. Plot expression without changing the saved default assay.
if (!(settings$lineage_assay %in% Assays(obj)))
  stop_bad("Missing lineage plotting assay: ", settings$lineage_assay)
previous_assay <- DefaultAssay(obj)
DefaultAssay(obj) <- settings$lineage_assay
lineage_markers <- settings$lineage_markers
missing <- setdiff(lineage_markers, rownames(obj[[settings$lineage_assay]]))
if (length(missing)) warning("Lineage markers absent from ", settings$lineage_assay,
  ": ", paste(missing, collapse = ", "), call. = FALSE)
lineage_markers <- setdiff(lineage_markers, missing)
if (length(lineage_markers)) {
  p_lineage <- FeatureDimPlot(obj, reduction = reduction, features = lineage_markers,
    slot = "data", pt.size = settings$plot_point_size, theme_use = "theme_blank",
    compare_features = TRUE, color_blend_mode = "blend",
    label = TRUE, label_insitu = TRUE,
    xlab = "wnnUMAP_1", ylab = "wnnUMAP_2", seed = 123)
} else {
  # Produce a readable report even when the reference uses different gene IDs.
  p_lineage <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = "No requested lineage markers found in the plotting assay") +
    ggplot2::theme_void()
}
ggplot2::ggsave(file.path(plot_dir, "Cell_lineage_umap.pdf"),
  plot = p_lineage, width = 5, height = 3.5)
DefaultAssay(obj) <- previous_assay
saveRDS(obj, outfile)
