#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) == 1L)
  obj <- readRDS(args[1])
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  Idents(obj) <- "wsnn_clusters"
  markers <- if (length(levels(Idents(obj))) > 1) {
    FindAllMarkers(obj, only.pos = TRUE, assay = "RNA", min.pct = 0.1,
                   logfc.threshold = 0.25, test.use = "wilcox")
  } else data.frame()
  # Apply significance to adjusted p-values explicitly.
  markers <- filter_significant_markers(markers)
  write_tsv(markers, "cluster_markers.tsv")
