#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) >= 4L)
  outfile <- args[1]; lo <- as.integer(args[2]); hi <- as.integer(args[3])
  gr <- GenomicRanges::reduce(do.call(c, lapply(args[-(1:3)],
                              function(f) rtracklayer::import(f, format = "BED"))),
                              ignore.strand = TRUE)
  gr <- gr[GenomicRanges::width(gr) > lo & GenomicRanges::width(gr) < hi]
  if (!length(gr)) stop_bad("No consensus peaks after width filtering")
  rtracklayer::export(gr, outfile, format = "BED")
