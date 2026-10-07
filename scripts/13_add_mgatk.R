#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  # Positional arguments are WNN, output RDS, output TSV and mgatk tables.
  # An optional final --copy-number-files section supplies step 10 tables.
  copy_marker <- which(args == "--copy-number-files")
  copy_files <- character()
  if (length(copy_marker)) {
    if (length(copy_marker) != 1L || copy_marker <= 4L || copy_marker == length(args))
      stop_bad("Use --copy-number-files after the mgatk tables, followed by copy-number TSVs")
    copy_files <- args[seq.int(copy_marker + 1L, length(args))]
    args <- args[seq_len(copy_marker - 1L)]
  }
  if (length(args) < 4L) stop_bad("Expected WNN RDS, output RDS, output TSV and mgatk tables")
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
    # TSV type inference turns an entirely missing numeric column into logical.
    # Retain numeric coverage/fraction metadata even when no cells cover a site.
    numeric_columns <- grep("^(mtDNA_(depth|count)$|Depth_|Informative_depth_|Allele_|[ACGT]_[0-9]|Heteroplasmy_)",
                            names(tab), value = TRUE)
    for (key in numeric_columns) {
      if (is.logical(tab[[key]]) && all(is.na(tab[[key]]))) tab[[key]] <- as.numeric(tab[[key]])
      values <- tab[[key]]
      if (!is.numeric(values) || any(!is.na(values) & (!is.finite(values) | values < 0)) ||
          (startsWith(key, "Heteroplasmy_") && any(values > 1, na.rm = TRUE)))
        stop_bad("Invalid mgatk measurement ", key, ": ", f)
    }
    tab
  })
  if (!all(vapply(tables, function(tab) identical(names(tab), names(tables[[1]])), logical(1))))
    stop_bad("mgatk tables have different columns; regenerate all samples with the current 09 script")
  tab <- do.call(rbind, tables)
  if (anyDuplicated(tab$cell_id)) stop_bad("Repeated mgatk cell IDs")
  if (!setequal(tab$cell_id, colnames(obj))) stop_bad("mgatk/Seurat cell ID mismatch")
  tab <- tab[match(colnames(obj), tab$cell_id), , drop = FALSE]
  if (!identical(tab$cell_id, colnames(obj))) stop_bad("mgatk cell order mismatch")
  for (key in c("sample_id", "barcode"))
    if (!identical(as.character(tab[[key]]), as.character(obj[[]][[key]])))
      stop_bad("mgatk/Seurat metadata mismatch: ", key)
  if (length(copy_files)) {
    # Copy-number tables use sample_id + original barcode, rather than Seurat's
    # stage-prefixed cell IDs. Never join on barcode alone: libraries reuse IDs.
    cell_key <- function(x) {
      identifiers <- x[, c("sample_id", "barcode"), drop = FALSE]
      if (anyNA(identifiers) || any(!nzchar(as.matrix(identifiers))) ||
          any(grepl("[[:cntrl:]]", as.matrix(identifiers))))
        stop_bad("Missing or invalid copy-number sample/barcode identifiers")
      paste(x$sample_id, x$barcode, sep = "\t")
    }
    target_keys <- cell_key(tab)
    if (anyDuplicated(target_keys)) stop_bad("Repeated sample/barcode pairs in final object")
    copy_tables <- lapply(copy_files, function(f) {
      values <- read.delim(f, check.names = FALSE,
        colClasses = c(sample_id = "character", barcode = "character"))
      estimates <- intersect(c("mtDNA_CN_normalised", "mtDNA_CN_normalised_enriched"), names(values))
      if (!all(c("sample_id", "barcode") %in% names(values)) ||
          !length(estimates) || anyDuplicated(names(values)) || !nrow(values))
        stop_bad("Missing identifiers/estimates or invalid copy-number columns: ", f)
      numeric_columns <- grep("^(mtDNA_CN_normalised|mtDNA_copy_estimate|nuclear_|mitochondrial_mean_depth|assumed_nuclear_ploidy)",
                              names(values), value = TRUE)
      for (key in numeric_columns) {
        # read.delim infers logical for an entirely NA column; preserve numeric
        # metadata even when all cells fail that source's coverage checks.
        if (is.logical(values[[key]]) && all(is.na(values[[key]]))) values[[key]] <- as.numeric(values[[key]])
        if (!is.numeric(values[[key]]) || any(!is.na(values[[key]]) &
            (!is.finite(values[[key]]) | values[[key]] < 0)))
          stop_bad("Invalid copy-number estimate ", key, ": ", f)
      }
      values
    })
    copy_keys <- unlist(lapply(copy_tables, cell_key), use.names = FALSE)
    if (anyDuplicated(copy_keys)) stop_bad("Repeated copy-number sample/barcode pairs")
    if (any(!copy_keys %in% target_keys)) stop_bad("Copy-number cells do not belong to the final object")
    # A source can be omitted for some samples; union columns and fill with NA.
    copy_columns <- unique(unlist(lapply(copy_tables, names), use.names = FALSE))
    copy_tables <- lapply(copy_tables, function(values) {
      for (key in setdiff(copy_columns, names(values))) values[[key]] <- NA
      values[, copy_columns, drop = FALSE]
    })
    copies <- do.call(rbind, copy_tables)
    copies <- copies[match(target_keys, copy_keys), , drop = FALSE]
    value_columns <- setdiff(copy_columns, c("sample_id", "barcode"))
    # Keep copy-number QC separate from any existing general cell QC metadata.
    renamed <- sub("^qc_status", "mtDNA_CN_qc_status", value_columns)
    if (anyDuplicated(renamed) || any(renamed %in% names(tab)))
      stop_bad("Copy-number columns conflict with mgatk metric columns")
    for (i in seq_along(value_columns)) tab[[renamed[i]]] <- copies[[value_columns[i]]]
    tab$mtDNA_CN_matched <- target_keys %in% copy_keys
  }
  rownames(tab) <- tab$cell_id
  obj <- AddMetaData(obj, metadata = tab[setdiff(names(tab),
    c("cell_id", "sample_id", "barcode"))])
  write_tsv(tab, tabfile)
  saveRDS(obj, outfile)
