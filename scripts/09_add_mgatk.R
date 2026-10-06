#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) >= 4L)
  obj <- readRDS(args[1]); outfile <- args[2]; tabfile <- args[3]
  tables <- lapply(args[-(1:3)], function(f) {
    # Preserve numeric-looking labels such as sample_id = "001" literally.
    tab <- read.delim(f, check.names = FALSE,
      colClasses = c(cell_id = "character", sample_id = "character", barcode = "character"))
    required <- c("cell_id", "sample_id", "barcode", "mtDNA_depth", "mtDNA_count",
                  "Heteroplasmy_5019", "Heteroplasmy_all_5019", "Depth_all_5019")
    if (!all(required %in% names(tab)) || anyDuplicated(names(tab)))
      stop_bad("Missing or duplicated mgatk metric columns: ", f)
    if (!nrow(tab) || anyNA(tab[, c("cell_id", "sample_id", "barcode")]) ||
        any(!nzchar(tab$cell_id)) || any(!nzchar(tab$sample_id)) || any(!nzchar(tab$barcode)))
      stop_bad("Empty mgatk table or missing identifiers: ", f)
    tab
  })
  if (!all(vapply(tables, function(tab) identical(names(tab), names(tables[[1]])), logical(1))))
    stop_bad("mgatk tables have different columns; regenerate all samples with the current 06c script")
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
