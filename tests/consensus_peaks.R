# Test actual union/BED coordinates and strict width boundaries without Signac.
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
run_checks <- function() {
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- file.path(root, c('a.bed', 'b.bed'))
  writeLines(c('chr1\t0\t20', 'chr1\t40\t70', 'chr1\t120\t10120'), files[1])
  writeLines(c('chr1\t70\t100', 'chr1\t10200\t10221', 'chr2\t0\t40'), files[2])
  out <- file.path(root, 'common.bed')
  expressions <- parse('scripts/02_consensus_peaks.R')[-(1:2)]
  run <- function(lo = '20', hi = '10000') {
    env$args <- c(out, lo, hi, files)
    for (expr in expressions) eval(expr, env)
  }
  run()
  tab <- read.table(out)[, 1:3]
  stopifnot(identical(tab[[1]], c('chr1', 'chr1', 'chr2')),
            identical(tab[[2]], c(40L, 10200L, 0L)),
            identical(tab[[3]], c(100L, 10221L, 40L)))
  run('20.0', '10000.0')
  for (pair in list(c('20.5', '10000'), c('20', '10000.5'), c('-1', '10000'),
                    c('20', '20'), c('NA', '10000'), c('20', '2147483648')))
    stopifnot(inherits(tryCatch(run(pair[1], pair[2]), error = identity), 'error'))
}
run_checks()
cat('All consensus-peak union/width checks passed.\n')
