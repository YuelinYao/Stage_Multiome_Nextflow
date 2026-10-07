#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
sample <- args[1]; arc <- args[2]; bed <- args[3]
message("03_qc_metrics.R arguments:\n",
        "  sample: ", sample, "\n",
        "  arc: ", arc, "\n",
        "  bed: ", bed, "\n",
        "  output TSV: ", args[4])
tab <- qc_table(sample, make_assays(arc, bed))
write_tsv(tab, args[4])
