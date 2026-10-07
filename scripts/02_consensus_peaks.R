#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) >= 4L)
  outfile <- args[1]; lo <- as.integer(args[2]); hi <- as.integer(args[3])
  message("02_consensus_peaks.R arguments:\n",
          "  outfile: ", outfile, "\n",
          "  peak_min: ", args[2], "\n",
          "  peak_max: ", args[3], "\n",
          "  input BED files:\n    ", paste(args[-(1:3)], collapse = "\n    "))
  if (anyNA(c(lo, hi)) || lo < 0L || hi <= lo)
    stop_bad("Peak widths require 0 <= min < max")
  gr <- GenomicRanges::reduce(do.call(c, lapply(args[-(1:3)],
                              function(f) rtracklayer::import(f, format = "BED"))),
                              ignore.strand = TRUE)
  gr <- gr[GenomicRanges::width(gr) > lo & GenomicRanges::width(gr) < hi]
  if (!length(gr)) stop_bad("No consensus peaks after width filtering")
  rtracklayer::export(gr, outfile, format = "BED")
