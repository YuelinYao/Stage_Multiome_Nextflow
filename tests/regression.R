# Run from the repository root; no Signac or real sequencing inputs required.
for (f in list.files("scripts", pattern = "[.]R$", full.names = TRUE)) parse(f)
env <- new.env(parent = globalenv())
# Load pure helpers without loading the unavailable scientific dependencies.
for (expr in parse("scripts/00_common.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-"))) eval(expr, env)
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
cfg <- jsonlite::fromJSON("examples/pipeline_settings.example.json", simplifyVector = FALSE)
config_file <- tempfile(fileext = ".json")
write_config <- function(x) {
  jsonlite::write_json(x, config_file, auto_unbox = TRUE, null = "null")
  config_file
}
valid <- env$read_pipeline_config(write_config(cfg))
stopifnot(identical(valid$review_QC_cluster$RNA_dims, 1:20))
stopifnot(identical(valid$ATAC_processing$min_cutoff, "q5"))
for (cutoff in list("q0", "q5", "q95.5", 0, 10)) {
  custom_cutoff <- cfg; custom_cutoff$ATAC_processing$min_cutoff <- cutoff
  stopifnot(isTRUE(all.equal(env$read_pipeline_config(write_config(custom_cutoff))$ATAC_processing$min_cutoff, cutoff)))
}
custom_cutoff$ATAC_processing$min_cutoff <- NULL
stopifnot(identical(env$read_pipeline_config(write_config(custom_cutoff))$ATAC_processing$min_cutoff, "q5"))
for (cutoff in list("q100", "q101", "q-5", "5", "invalid", -1, 1.5, TRUE)) {
  custom_cutoff$ATAC_processing$min_cutoff <- cutoff
  expect_error(env$read_pipeline_config(write_config(custom_cutoff)))
}
# Every section can differ; none should inherit another step's settings.
independent <- cfg
for (section in c("review_QC_cluster", "doublet", "RNA_processing", "ATAC_processing", "WNN")) {
  independent[[section]]$algorithm <- match(section, c("review_QC_cluster", "doublet", "RNA_processing", "ATAC_processing"), nomatch = 1)
  independent[[section]]$resolution <- match(section, names(independent)) / 10
}
independent$review_QC_cluster$RNA_dims <- "1:12"
independent$doublet$PC_dims <- "1:15"
independent$RNA_processing$RNA_dims <- "1:25"
independent$WNN$RNA_dims <- "1:35"
independent$WNN$ATAC_dims <- "2:45"
independent$ATAC_processing$ATAC_integrate_dims <- "1:45"
separate <- env$read_pipeline_config(write_config(independent))
stopifnot(identical(separate$review_QC_cluster$RNA_dims, 1:12),
          identical(separate$doublet$PC_dims, 1:15),
          identical(separate$RNA_processing$RNA_dims, 1:25),
          identical(separate$WNN$RNA_dims, 1:35),
          separate$doublet$algorithm == 2,
          separate$RNA_processing$algorithm == 3,
          separate$ATAC_processing$algorithm == 4,
          separate$WNN$algorithm == 1)
bad <- independent; bad$ATAC_processing$ATAC_integrate_dims <- "1:40"
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$ATAC_processing$ATAC_integrate_dims <- "2:40"
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$doublet$resolution <- 0
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$clustering <- list(algorithm = 3)
expect_error(env$read_pipeline_config(write_config(bad)))
stopifnot(valid$peaks$peak_min == 20, valid$peaks$peak_max == 10000)
legacy <- cfg; legacy$peaks <- NULL
stopifnot(env$read_pipeline_config(write_config(legacy))$peaks$peak_min == 20)
bad <- cfg; bad$peaks$peak_min <- 10000
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$peaks$peak_max <- 12.5
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$sample_overrides <- list(list(nCount_RNA = list(min = 1, max = 2)))
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$thresholds$nCount_RNA <- list(min = 2, max = 1)
expect_error(env$read_pipeline_config(write_config(bad)))
bad <- cfg; bad$review_QC_cluster$RNA_dims <- "20:1"
expect_error(env$read_pipeline_config(write_config(bad)))
tab <- data.frame(barcode = c("lower", "upper", "below", "missing", "infinite"))
limits <- setNames(lapply(env$qc_metric_names, function(x) list(min = 1, max = 2)),
                   env$qc_metric_names)
for (metric in env$qc_metric_names) tab[[metric]] <- c(1, 2, 0, NA, Inf)
stopifnot(identical(env$apply_qc(tab, limits), c("lower", "upper")))

# Execute the mgatk extractor with a mock ReadMGATK and a plain input object.
obj_file <- tempfile(); out_file <- tempfile()
obj <- matrix(0, 1, 2, dimnames = list("gene", c("sample_present", "sample_missing")))
class(obj) <- c("mock_cells", class(obj))
`$.mock_cells` <- function(x, name) c("present", "missing")
saveRDS(obj, obj_file)
env$ReadMGATK <- function(dir) {
  positions <- c(5019, 15200, 16232)
  rows <- unlist(lapply(positions, function(pos)
    paste(rep(c("A", "G", "C", "T"), each = 2), pos,
          rep(c("fwd", "rev"), 4), sep = "-")))
  counts <- matrix(0, nrow = length(rows), ncol = 1,
                   dimnames = list(rows, "present"))
  counts[paste(rep(c("A", "G"), each = 2), 5019,
               rep(c("fwd", "rev"), 2), sep = "-"), 1] <- c(3, 3, 2, 2)
  list(refallele = data.frame(pos = positions,
    ref = c("A", "A", "A")),
    counts = counts, depth = data.frame(depth = 100, row.names = "present"))
}
env$args <- c("sample", obj_file, "mock", out_file)
exprs <- parse("scripts/06c_mgatk_metrics.R")
for (expr in exprs[-(1:2)]) eval(expr, env)
result <- read.delim(out_file)
stopifnot(result$Heteroplasmy_5019[1] == 0.4, result$Depth_5019[1] == 10,
          result$mtDNA_depth[1] == 100, is.na(result$Heteroplasmy_5019[2]),
          is.na(result$mtDNA_depth[2]))
unlink(c(config_file, obj_file, out_file))
stopifnot(env$nextgem_multiome_doublet_rate(5000) == 0.039,
          abs(env$nextgem_multiome_doublet_rate(5500) - 0.0425) < 1e-12)
expect_error(env$nextgem_multiome_doublet_rate(10001))
expect_error(env$nextgem_multiome_doublet_rate(499))
stopifnot(identical(valid$doublet$expected_rate, "auto"))
stopifnot(identical(valid$doublet$PC_dims, 1:20))
custom <- cfg; custom$doublet$PC_dims <- "1:30"
stopifnot(identical(env$read_pipeline_config(write_config(custom))$doublet$PC_dims, 1:30))
custom$doublet$PC_dims <- NULL
stopifnot(identical(env$read_pipeline_config(write_config(custom))$doublet$PC_dims, 1:20))
for (dims in c("2:20", "1:51", "20:1", "invalid")) {
  custom$doublet$PC_dims <- dims
  expect_error(env$read_pipeline_config(write_config(custom)))
}
override <- cfg; override$doublet$sample_rates <- list(WT_1 = 0.04)
stopifnot(env$read_pipeline_config(write_config(override), "WT_1")$doublet$expected_rate == 0.04)
unlink(config_file)
cat("All regression checks passed.\n")
