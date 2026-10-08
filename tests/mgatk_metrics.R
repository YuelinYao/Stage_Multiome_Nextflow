# Synthetic sparse mgatk counts exercise actual 09 and actual Seurat metadata
# attachment in 13. ReadMGATK is mocked; production mgatk file parsing is not.
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
  positions <- c(1781, 1866, 3009, 3823, 5019, 5024, 13614, 13715, 15200, 16232)
  refs <- c('C', 'A', 'G', 'T', 'A', 'C', 'C', 'C', 'A', 'A')
  alts <- c('T', 'G', 'T', 'C', 'G', 'T', 'T', 'T', 'G', 'T')
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
  expressions <- parse('scripts/09_mgatk_metrics.R')[-(1:2)]
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
    tab$mtDNA_depth[1] == 30, tab$mtDNA_count[1] == 120,
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
    tab$Heteroplasmy_5024[1] == 0.7,
    tab$C_5024[1] == 3, tab$T_5024[1] == 7,
    tab$Depth_all_5024[1] == 84,
    tab$Informative_depth_all_5024[1] == 70,
    tab$Allele_ref_all_5024[1] == 41, tab$Allele_mut_all_5024[1] == 29,
    abs(tab$Heteroplasmy_all_5024[1] - 29/70) < 1e-12,
    is.na(tab$Heteroplasmy_all_5024_covered[4]),
    !any(grepl('all_5024_reciprocal', names(tab))))
  # Every per-site fraction follows its stated allele definition.
  for (pos in positions) stopifnot(tab[[paste0('Heteroplasmy_', pos)]][1] == 0.7)

  # Unequal coverage verifies that pools sum counts, not site fractions.
  original_counts <- mgatk$counts
  mgatk$counts[, 'covered'] <- 0
  mgatk$counts['A-5019-fwd', 'covered'] <- 10
  mgatk$counts['G-15200-rev', 'covered'] <- 90
  # A reciprocal 5024 haplotype has reference bases at five linked sites,
  # but alternate bases at 5024 and 13715, exactly as in YY_06.
  pool_positions <- c(1781, 1866, 3009, 3823, 5024, 13614, 13715)
  for (pos in pool_positions) {
    i <- match(pos, positions)
    base <- if (pos %in% c(5024, 13715)) alts[i] else refs[i]
    mgatk$counts[paste(base, pos, 'fwd', sep = '-'), 'covered'] <- 10
  }
  weighted <- run()
  stopifnot(weighted$Heteroplasmy_all_5019[1] == 0.9,
    weighted$Heteroplasmy_all_5024[1] == 1,
    weighted$Heteroplasmy_all_5024_covered[1] == 1)
  # The opposite haplotype must give zero, whatever the sample label.
  mgatk$counts[, 'covered'] <- 0
  for (pos in pool_positions) {
    i <- match(pos, positions)
    base <- if (pos %in% c(5024, 13715)) refs[i] else alts[i]
    mgatk$counts[paste(base, pos, 'fwd', sep = '-'), 'covered'] <- 10
  }
  for (label in c('5019_example', '5024_example', 'WT_example')) {
    env$args <- c(label, input, 'mock', output)
    for (expr in expressions) eval(expr, env)
    result <- read.delim(output)
    stopifnot(result$Heteroplasmy_all_5024[1] == 0,
              result$Heteroplasmy_all_5024_covered[1] == 0)
  }
  mgatk$counts <- original_counts

  # Actual 13 must attach the expanded table by cell ID, not input row order.
  tab <- tab[c(4, 2, 1, 3), ]; env$write_tsv(tab, output)
  final <- file.path(task_dir, 'wnn_mgatk.rds'); combined <- file.path(task_dir, 'combined.tsv')
  env$args <- c(input, final, combined, output)
  for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env)
  attached <- readRDS(final)
  stopifnot(identical(colnames(attached), colnames(obj)),
    identical(as.character(attached$sample_id), rep('001', ncol(obj))),
    attached$Heteroplasmy_all_5019[1] == 0.7,
    attached$Heteroplasmy_5024[1] == 0.7,
    abs(attached$Heteroplasmy_all_5024[1] - 29/70) < 1e-12, is.na(attached$Heteroplasmy_all_5019[2]))
  # All-NA numeric fractions retain numeric metadata after a TSV round trip.
  all_missing <- tab
  all_missing$Heteroplasmy_15200_covered <- NA_real_
  env$write_tsv(all_missing, output)
  for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env)
  missing_result <- readRDS(final)
  stopifnot(is.numeric(missing_result$Heteroplasmy_15200_covered),
            all(is.na(missing_result$Heteroplasmy_15200_covered)))
  env$write_tsv(tab, output)
  attachment_args <- env$args
  # Step 13 imports both copy-number sources without disturbing heteroplasmy.
  copy_file <- file.path(task_dir, '001.mtDNA_copy_number.tsv')
  copies <- data.frame(sample_id = '001', barcode = c('low', 'covered', 'other_only'),
    mtDNA_CN_normalised = c(NA, 200, 60),
    mtDNA_CN_normalised_enriched = c(20, 800, NA),
    qc_status = c('low_mitochondrial_coverage', 'pass', 'pass'),
    qc_status_enriched = c('pass', 'pass', 'low_mitochondrial_coverage'))
  env$write_tsv(copies, copy_file)
  attach_copies <- function(files = copy_file) {
    env$args <- c(attachment_args, '--copy-number-files', files)
    for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env)
    readRDS(final)
  }
  annotated <- attach_copies()
  stopifnot(identical(as.numeric(annotated$mtDNA_CN_normalised), c(200, NA_real_, 60, NA_real_)),
    identical(as.numeric(annotated$mtDNA_CN_normalised_enriched), c(800, NA_real_, NA_real_, 20)),
    identical(unname(annotated$mtDNA_CN_matched), c(TRUE, FALSE, TRUE, TRUE)),
    annotated$mtDNA_CN_qc_status[4] == 'low_mitochondrial_coverage',
    annotated$Heteroplasmy_all_5019[1] == 0.7,
    identical(as.numeric(read.delim(combined)$mtDNA_CN_normalised), as.numeric(annotated$mtDNA_CN_normalised)))
  # All-missing estimates must still be numeric metadata.
  all_missing <- copies; all_missing$mtDNA_CN_normalised_enriched <- NA_real_
  env$write_tsv(all_missing, copy_file)
  annotated <- attach_copies()
  stopifnot(is.numeric(annotated$mtDNA_CN_normalised_enriched),
    all(is.na(annotated$mtDNA_CN_normalised_enriched)))
  env$write_tsv(copies, copy_file)
  expect_error(attach_copies(c(copy_file, copy_file)))
  bad_copies <- copies; bad_copies$sample_id <- 'different_sample'
  env$write_tsv(bad_copies, copy_file); expect_error(attach_copies())
  bad_copies <- copies; bad_copies$mtDNA_CN_normalised[1] <- -1
  env$write_tsv(bad_copies, copy_file); expect_error(attach_copies())
  bad_copies <- copies; bad_copies$barcode[1] <- NA_character_
  env$write_tsv(bad_copies, copy_file); expect_error(attach_copies())
  # Reused barcodes across samples must resolve using BOTH identifiers.
  second_obj <- obj
  colnames(second_obj) <- paste0('second_', bcs)
  second_obj$sample_id <- '002'
  stage_obj <- merge(obj, second_obj)
  stage_input <- file.path(task_dir, 'stage.rds'); saveRDS(stage_obj, stage_input)
  second_tab <- tab; second_tab$cell_id <- paste0('second_', second_tab$barcode)
  second_tab$sample_id <- '002'
  second_metrics <- file.path(task_dir, 'second_metrics.tsv')
  env$write_tsv(second_tab, second_metrics)
  second_copies <- data.frame(sample_id = '002', barcode = 'covered', mtDNA_CN_normalised = 123)
  second_file <- file.path(task_dir, '002.mtDNA_copy_number.tsv')
  env$write_tsv(second_copies, second_file); env$write_tsv(copies, copy_file)
  env$args <- c(stage_input, final, combined, output, second_metrics,
    '--copy-number-files', second_file, copy_file)
  for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env)
  annotated <- readRDS(final)
  stopifnot(annotated[[]]['alias_covered', 'mtDNA_CN_normalised'] == 200,
    annotated[[]]['second_covered', 'mtDNA_CN_normalised'] == 123,
    is.na(annotated[[]]['second_covered', 'mtDNA_CN_normalised_enriched']))
  env$args <- attachment_args
  # Malformed tables must fail clearly instead of attaching partial metadata.
  bad <- tab; bad$cell_id[1] <- NA_character_; env$write_tsv(bad, output)
  expect_error(for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env))
  bad <- tab; bad$Heteroplasmy_all_5019 <- NULL; env$write_tsv(bad, output)
  env$args <- attachment_args
  expect_error(for (expr in parse('scripts/13_add_mgatk.R')[-(1:2)]) eval(expr, env))

  original <- mgatk
  mgatk$refallele$ref[mgatk$refallele$pos == 15200] <- 'C'
  expect_error(run()); mgatk <- original
  mgatk$refallele$ref[mgatk$refallele$pos == 5024] <- 'A'
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
