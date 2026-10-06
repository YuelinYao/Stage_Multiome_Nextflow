# Synthetic checks for the actual step 04 script and shared plotting helper.
# Run from the repository root: Rscript tests/qc_review.R
run_checks <- function() {
  original_dir <- getwd()
  script <- normalizePath('scripts/04_qc_review.R')
  supplied_json <- normalizePath('examples/pipeline_settings.example.json')
  env <- new.env(parent = globalenv())
  for (expr in parse('scripts/00_common.R'))
    if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
  root <- tempfile('qc-review-'); dir.create(root)
  on.exit({ setwd(original_dir); unlink(root, recursive = TRUE) }, add = TRUE)
  setwd(root)
  fixture <- function(id, values) {
    tab <- data.frame(sample_id = id, barcode = sprintf('%03d', seq_along(values)))
    for (metric in env$qc_metric_names) tab[[metric]] <- values
    tab
  }
  a <- fixture('A', c(1, 2, 3, 4, NA_real_, Inf))
  b <- fixture('B', c(10, 20, 30, 40, 50, 60))
  env$write_tsv(a, 'a.tsv'); env$write_tsv(b, 'b.tsv')
  file.copy(supplied_json, 'pipeline_settings.json')
  original_json <- readBin('pipeline_settings.json', 'raw', file.info('pipeline_settings.json')$size)
  env$args <- c('summary.tsv', 'mad.tsv', 'stage.pdf', 'b.tsv', 'a.tsv')
  # Skip only the source() and commandArgs() statements; execute the real code.
  expressions <- parse(script)[-(1:2)]
  for (expr in expressions) eval(expr, env)
  combined <- read.delim('all_cells_qc_metrics.tsv', colClasses = c(barcode='character'))
  stopifnot(nrow(combined) == 12, combined$sample_id[1] == 'A', combined$barcode[1] == '001')
  summary <- read.delim('summary.tsv')
  selected <- subset(summary, sample_id == 'A' & metric == 'nCount_RNA')
  stopifnot(nrow(summary) == 24, selected$n == 4, selected$median == 2.5,
    selected$mad_scaled == 1.4826,
    isTRUE(all.equal(selected$mad_3_lower, 2.5 - 3 * 1.4826)))
  pooled <- subset(summary, sample_id == 'ALL_STAGE' & metric == 'nCount_RNA')
  stopifnot(pooled$n == 10, pooled$median == median(c(1:4, seq(10,60,10))))
  candidates <- read.delim('mad.tsv')
  stopifnot(nrow(candidates) == 72, setequal(candidates$mad_multiplier, c(3,5,7)))
  stopifnot(file.info('stage.pdf')$size > 0,
    identical(list.files(pattern = '[.]json$'), 'pipeline_settings.json'),
    identical(readBin('pipeline_settings.json', 'raw', file.info('pipeline_settings.json')$size),
              original_json))
  # Malformed and duplicated input data should fail before publishing results.
  bad <- a; bad$nCount_RNA <- 'invalid'; env$write_tsv(bad, 'bad.tsv')
  env$args <- c(env$args[1:3], 'bad.tsv')
  stopifnot(inherits(tryCatch({for(expr in expressions) eval(expr, env)}, error=identity), 'error'))
  env$args <- c(env$args[1:3], 'a.tsv', 'a.tsv')
  stopifnot(inherits(tryCatch({for(expr in expressions) eval(expr, env)}, error=identity), 'error'))
}
run_checks()
cat('All QC review checks passed.\n')
