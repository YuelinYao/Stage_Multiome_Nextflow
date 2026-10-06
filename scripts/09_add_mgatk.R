#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) >= 4L)
  obj <- readRDS(args[1]); outfile <- args[2]; tabfile <- args[3]
  tables <- lapply(args[-(1:3)], function(f) read.delim(f, check.names = FALSE))
  tab <- do.call(rbind, tables)
  if (anyDuplicated(tab$cell_id)) stop_bad("Repeated mgatk cell IDs")
  if (!setequal(tab$cell_id, colnames(obj))) stop_bad("mgatk/Seurat cell ID mismatch")
  tab <- tab[match(colnames(obj), tab$cell_id), , drop = FALSE]
  if (!identical(tab$cell_id, colnames(obj))) stop_bad("mgatk cell order mismatch")
  for (key in c("sample_id", "barcode"))
    if (!identical(as.character(tab[[key]]), as.character(obj[[]][[key]])))
      stop_bad("mgatk/Seurat metadata mismatch: ", key)
  rownames(tab) <- tab$cell_id
  obj <- AddMetaData(obj, metadata = tab[setdiff(names(tab),
    c("cell_id", "sample_id", "barcode"))])
  write_tsv(tab, tabfile)
  saveRDS(obj, outfile)
