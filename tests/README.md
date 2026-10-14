# Tests

Synthetic checks of the pipeline scripts, settings and Nextflow caching. They use small synthetic data and stand-ins for slow or unavailable tools (for example DoubletFinder, `ReadMGATK` or cellXY), so they check the pipeline's own logic and data handling, not the biology or the external packages. See the [Tests section of the main README](../README.md#tests) for how to run them.

| Test | Checks |
| --- | --- |
| `consensus_peaks.R` | Step 02: peak union, BED coordinates and width limits |
| `qc_review.R` | Step 04: combined QC tables, MAD bounds and plots |
| `reviewed_clusters.R` | Steps 05 → 06: saved cluster labels, QC-passing cells and stale-table detection |
| `doubletfinder_stage.R` | Step 06: QC, cluster removal, DoubletFinder routing, cellXY male-female doublet removal (`remove_mf_doublets` on and off) and the report |
| `rna_integration.R` | Step 07: per-sample processing, cell-cycle scoring, SCT and CCA integration, annotation table |
| `atac_integration.R` | Step 08: LSI, reciprocal-LSI integration, feature assays carried into both objects, motif-ID-to-TF-name renaming |
| `mgatk_metrics.R` | Steps 09 and 13: heteroplasmy, depth and pooled metrics, chrM coverage table and plot with enriched, ordinary or both mgatk runs, variant flags, tables and plot with default and custom `MGATK` thresholds, attachment to the final object |
| `mtdna_copy_number.R` | Step 10: reference selection, base counting, barcode translation, both mgatk sources and QC statuses |
| `copy_number_batch.R` | Step 10: batched base counter matches the original implementation (`--benchmark` for timing) |
| `wnn_settings.R` | Step 12: WNN settings, reduction choice and copying of step 08 assays |
| `modalities_markers.R` | Step 14 markers (combined and per-cluster files) and the legacy step 11 helper |
| `regression.R` | Shared helpers: independent section settings, peak counting across Signac versions, mgatk extraction |
| `process_settings.R` | R settings reader on the task JSONs (including `MGATK`), defaults and invalid values |
| `pipeline_settings.groovy` | Groovy settings validation, defaults and per-task settings |
| `settings_cache.py` | Real Nextflow on mock data: which tasks rerun after each kind of edit, sample-ID checks, step overlap; also validates every generated task JSON with the R reader |
| `script_hash.py` | Real Nextflow: editing an R script reruns tasks |
