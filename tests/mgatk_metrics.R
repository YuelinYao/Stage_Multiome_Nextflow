# Synthetic sparse mgatk counts exercise actual 06c and actual Seurat metadata
# attachment in 09. ReadMGATK is mocked; production mgatk file parsing is not.
suppressPackageStartupMessages({ library(Matrix); library(Seurat) })
env <- new.env(parent = globalenv())
for (expr in parse('scripts/00_common.R'))
  if (is.call(expr) && identical(expr[[1]], as.name('<-'))) eval(expr, env)
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
}
run_checks <- function() {
  task_dir <- tempfile(); dir.create(task_dir)
  on.exit(unlink(task_dir, recursive = TRUE), add = TRUE)
  positions <- c(5019, 15200, 16232)
  refs <- c('A', 'A', 'A')
  alts <- c('G', 'G', 'T')
  rows <- unlist(lapply(positions, function(pos)
    paste(rep(c('A', 'C', 'G', 'T'), each = 2), pos,
          rep(c('fwd', 'rev'), 4), sep = '-')))
  counts <- matrix(0, length(rows), 3,
                   dimnames = list(rows, c('other_only', 'low', 'covered')))
  for (i in seq_along(positions)) {
    counts[paste(refs[i], positions[i], c('fwd', 'rev'), sep = '-'), 'covered'] <- c(1, 2)
    counts[paste(alts[i], positions[i], c('fwd', 'rev'), sep = '-'), 'covered'] <- c(3, 4)
    other <- setdiff(c('A', 'C', 'G', 'T'), c(refs[i], alts[i]))[1]
    counts[paste(other, positions[i], 'fwd', sep = '-'), 'covered'] <- 2
    # Low informative coverage cannot be rescued by many other-base calls.
    counts[paste(refs[i], positions[i], 'fwd', sep = '-'), 'low'] <- 1
    counts[paste(alts[i], positions[i], 'fwd', sep = '-'), 'low'] <- 1
    counts[paste(other, positions[i], 'fwd', sep = '-'), 'low'] <- 100
  }
  counts['C-5019-fwd', 'other_only'] <- 20
  # Scramble BOTH dimensions, independently of the Seurat object's order.
  set.seed(12)
  mgatk <- list(counts = as(counts[sample(rows), c('low', 'other_only', 'covered')], 'dgCMatrix'),
    refallele = data.frame(pos = positions, ref = refs),
    depth = data.frame(depth = c(10, 30, 20), row.names = c('low', 'covered', 'other_only')))
  bcs <- c('covered', 'missing', 'other_only', 'low')
  obj <- CreateSeuratObject(matrix(1, 2, 4,
    dimnames = list(c('gene1', 'gene2'), paste0('alias_', bcs))))
  obj$barcode <- bcs; obj$sample_id <- '001'
  input <- file.path(task_dir, 'cells.rds'); output <- file.path(task_dir, 'metrics.tsv')
  saveRDS(obj, input)
  env$ReadMGATK <- function(dir) mgatk
  expressions <- parse('scripts/06c_mgatk_metrics.R')[-(1:2)]
  run <- function() {
    env$args <- c('001', input, 'mock', output)
    for (expr in expressions) eval(expr, env)
    read.delim(output, colClasses = c(cell_id = "character", sample_id = "character", barcode = "character"))
  }
  tab <- run()
  stopifnot(identical(tab$sample_id, rep("001", nrow(tab))))
  yy06_names <- c('mtDNA_depth', 'mtDNA_count',
    paste0('Heteroplasmy_', positions), paste0('Depth_', positions),
    'Heteroplasmy_all_5019', 'Depth_all_5019')
  stopifnot(all(yy06_names %in% names(tab)),
    !any(grepl('^hf_|^depth_', names(tab))),
    tab$Allele_ref_5019[1] == 3, tab$Allele_mut_5019[1] == 7,
    tab$Allele_ref_all_5019[1] == 9, tab$Allele_mut_all_5019[1] == 21)
  stopifnot(identical(tab$cell_id, colnames(obj)), identical(tab$barcode, bcs),
    identical(tab$mgatk_matched, c(TRUE, FALSE, TRUE, TRUE)),
    tab$mtDNA_depth[1] == 30, tab$mtDNA_count[1] == 36,
    all(is.na(tab[2, setdiff(names(tab), c('cell_id', 'sample_id', 'barcode', 'mgatk_matched'))])),
    tab$Depth_5019[1] == 12, tab$Informative_depth_5019[1] == 10,
    tab$A_5019[1] == 3, tab$G_5019[1] == 7,
    tab$Heteroplasmy_5019[1] == 0.7, tab$Heteroplasmy_5019_covered[1] == 0.7,
    tab$Heteroplasmy_5019_fwd[1] == 0.75,
    abs(tab$Heteroplasmy_5019_rev[1] - 2/3) < 1e-12,
    tab$Depth_5019[3] == 20, is.na(tab$Heteroplasmy_5019[3]),
    tab$Heteroplasmy_5019[4] == 0.5, is.na(tab$Heteroplasmy_5019_covered[4]),
    tab$Depth_all_5019[1] == 36, tab$Informative_depth_all_5019[1] == 30,
    tab$Heteroplasmy_all_5019_covered[1] == 0.7,
    tab$Depth_all_5019[4] == 306, is.na(tab$Heteroplasmy_all_5019_covered[4]),
    !any(grepl('5024', names(tab))))
  # Every per-site fraction follows its stated allele definition.
  for (pos in positions) stopifnot(tab[[paste0('Heteroplasmy_', pos)]][1] == 0.7)

  # Actual 09 must attach the expanded table by cell ID, not input row order.
  tab <- tab[c(4, 2, 1, 3), ]; env$write_tsv(tab, output)
  final <- file.path(task_dir, 'wnn_mgatk.rds'); combined <- file.path(task_dir, 'combined.tsv')
  env$args <- c(input, final, combined, output)
  for (expr in parse('scripts/09_add_mgatk.R')[-(1:2)]) eval(expr, env)
  attached <- readRDS(final)
  stopifnot(identical(colnames(attached), colnames(obj)),
    identical(as.character(attached$sample_id), rep('001', ncol(obj))),
    attached$Heteroplasmy_all_5019[1] == 0.7, is.na(attached$Heteroplasmy_all_5019[2]))
  attachment_args <- env$args
  # Malformed tables must fail clearly instead of attaching partial metadata.
  bad <- tab; bad$cell_id[1] <- NA_character_; env$write_tsv(bad, output)
  expect_error(for (expr in parse('scripts/09_add_mgatk.R')[-(1:2)]) eval(expr, env))
  bad <- tab; bad$Heteroplasmy_all_5019 <- NULL; env$write_tsv(bad, output)
  env$args <- attachment_args
  expect_error(for (expr in parse('scripts/09_add_mgatk.R')[-(1:2)]) eval(expr, env))

  original <- mgatk
  mgatk$refallele$ref[mgatk$refallele$pos == 15200] <- 'C'
  expect_error(run()); mgatk <- original
  rownames(mgatk$counts)[2] <- rownames(mgatk$counts)[1]
  expect_error(run()); mgatk <- original
  colnames(mgatk$counts)[2] <- colnames(mgatk$counts)[1]
  expect_error(run()); mgatk <- original
  mgatk$counts <- mgatk$counts[setdiff(rownames(mgatk$counts), 'G-5019-fwd'), ]
  expect_error(run()); mgatk <- original
  mgatk$counts@x[1] <- -1
  expect_error(run()); mgatk <- original
  mgatk$depth$depth[1] <- -1
  expect_error(run()); mgatk <- original
  colnames(mgatk$counts) <- paste0('no_match_', seq_len(ncol(mgatk$counts)))
  expect_error(run())
}
run_checks()
cat('All mgatk multisite, pooled-allele and metadata checks passed.\n')
