nextflow.enable.dsl = 2

params.samples = null
params.stage = null
params.outdir = 'results'
params.macs3_path = null
params.json = null
params.cluster_review_dir = null
params.copy_number_dir = null
params.phase = 'review'
params.peak_min = null
params.peak_max = null

// Quote every shell argument, including paths containing spaces or apostrophes.
def shellQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
// Include R sources in task scripts so -resume notices helper/script edits.
def sourceDigest = java.security.MessageDigest.getInstance('SHA-256')
file("${projectDir}/scripts").toFile().listFiles()
  .findAll { it.name.endsWith('.R') }.sort { it.name }.each {
    sourceDigest.update(it.name.getBytes('UTF-8'))
    sourceDigest.update(it.bytes)
  }
def scriptsHash = sourceDigest.digest().encodeHex().toString()

if (!params.samples || !params.stage) error 'Provide --samples manifest.csv and --stage E8.5'
if (!(params.phase in ['review', 'cluster_review', 'analysis'])) error '--phase must be review, cluster_review or analysis'
if (!(params.stage.toString() ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Stage must be a safe directory name'
// Read only peak settings here: a review JSON can still have unset QC bounds.
def peakSettings = [:]
def pipelineSettings = [:]
if (params.json) {
  if (!file(params.json).exists()) error "Pipeline JSON not found: ${params.json}"
  def settings = new groovy.json.JsonSlurper().parse(file(params.json).toFile())
  if (!(settings instanceof Map)) error 'Pipeline JSON must be an object'
  def unknownSettings = settings.keySet() - ['peaks', 'note', 'thresholds', 'sample_overrides',
    'review_QC_cluster', 'doublet', 'RNA_processing', 'ATAC_processing', 'WNN', 'copy_number', 'clustering']
  if (unknownSettings) error "Unknown pipeline JSON sections: ${unknownSettings.join(', ')}"
  pipelineSettings = settings
  if (settings.peaks != null) {
    if (!(settings.peaks instanceof Map) ||
        !settings.peaks.keySet().every { it in ['peak_min', 'peak_max'] })
      error 'Pipeline JSON peaks must be an object containing peak_min and/or peak_max'
    peakSettings = settings.peaks
  }
}
// Explicit command-line values override JSON; omitted settings use defaults.
params.peak_min = params.peak_min != null ? params.peak_min : (peakSettings.peak_min != null ? peakSettings.peak_min : 20)
params.peak_max = params.peak_max != null ? params.peak_max : (peakSettings.peak_max != null ? peakSettings.peak_max : 10000)
if (!(params.peak_min.toString() ==~ /[0-9]+/) || !(params.peak_max.toString() ==~ /[0-9]+/) ||
    params.peak_min.toString().toLong() >= params.peak_max.toString().toLong())
  error 'Peak widths must be nonnegative integers with peak_min < peak_max'
if (!params.macs3_path) error 'Provide --macs3_path /absolute/path/to/macs3'
if (!java.nio.file.Paths.get(params.macs3_path.toString()).isAbsolute()) error '--macs3_path must be absolute'
if (!java.nio.file.Files.isRegularFile(file(params.macs3_path)) ||
    !java.nio.file.Files.isExecutable(file(params.macs3_path)))
  error "MACS3 executable not found or not executable: ${params.macs3_path}"
if (params.phase in ['cluster_review', 'analysis'] && (!params.json || !file(params.json).exists()))
  error 'Cluster review and analysis require --json /absolute/path/to/pipeline_settings.json'
if (params.enable_conda) {
  if (!params.conda_env) error 'The conda profile requires --conda_env /shared/path/to/env (or environment.yml)'
  if (!file(params.conda_env).exists()) error "Conda environment or YAML not found: ${params.conda_env}"
}

process PEAK_SAMPLE {
  tag "$id"
  cpus 4
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/peaks/sample", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk_enriched)
  output:
  tuple val(id), path("${id}.bed")
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/01_peaks.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(params.macs3_path)} '${id}.bed'
  """
}

process COMMON_PEAKS {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/peaks", mode: 'copy'
  input:
  path beds
  output:
  path 'common_peaks.bed'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/02_consensus_peaks.R")} common_peaks.bed ${params.peak_min} ${params.peak_max} ${beds.collect { shellQuote(it) }.join(' ')}
  """
}

process QC_METRICS {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/review", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk_enriched), path(peaks)
  output:
  tuple val(id), path("${id}.metrics.tsv"), path("${id}.distributions.pdf")
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/03_qc_metrics.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} '${id}.metrics.tsv' '${id}.distributions.pdf'
  """
}

process QC_REVIEW {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/qc/review", mode: 'copy'
  input:
  path metric_files
  output:
  path 'qc_summary.tsv'
  path 'qc_mad_candidates.tsv'
  path 'all_cells_qc_metrics.tsv'
  path 'stage_qc_distributions.pdf'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/04_qc_review.R")} qc_summary.tsv qc_mad_candidates.tsv stage_qc_distributions.pdf ${metric_files.collect { shellQuote(it) }.join(' ')}
  """
}

process PRE_DOUBLET_REVIEW {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/cluster_review", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk_enriched), path(peaks), path(json)
  output:
  path "${id}.cluster_cells.tsv"
  path "${id}.cluster_summary.tsv"
  path "${id}.cluster_qc.pdf"
  path "${id}.cluster_markers.tsv"
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/05_review_clusters.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} ${shellQuote(json)}
  """
}

process QC_FILTER_DOUBLET {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/filtered", mode: 'copy', pattern: '*.{qc.tsv,pdf}'
  input:
  tuple val(id), val(arc), val(mgatk_enriched), path(peaks), path(json), path(review)
  output:
  tuple val(id), val(mgatk_enriched), path("${id}.rna.qc.rds"), path("${id}.atac.qc.rds"), path("${id}.qc.tsv"), path("${id}_Doublet.pdf")
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06_qc_doubletfinder.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} ${shellQuote(json)} '${id}.rna.qc.rds' '${id}.atac.qc.rds' '${id}.qc.tsv' ${shellQuote(review)}
  """
}

process MERGE_INTEGRATE_RNA {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'rna_*'
  input:
  path objects
  path json
  output:
  path 'rna_cell_annotations.tsv', emit: annotations
  path 'rna_unintegrated.rds', emit: unintegrated
  path 'rna_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/07_preprocess_rna.R")} ${shellQuote(json)} ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MERGE_INTEGRATE_ATAC {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'atac_*.rds'
  input:
  path objects
  path json
  output:
  path 'atac_unintegrated.rds', emit: unintegrated
  path 'atac_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/08_preprocess_atac.R")} ${shellQuote(json)} ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MGATK_METRICS {
  tag "$id"
  cpus 2
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/mgatk", mode: 'copy'
  input:
  tuple val(id), val(mgatk_enriched), path(qc)
  output:
  tuple val(id), path("${id}.mtDNA.tsv")
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/09_mgatk_metrics.R")} ${shellQuote(id)} ${shellQuote(qc)} ${shellQuote(mgatk_enriched)} '${id}.mtDNA.tsv'
  """
}

// Step 10 follows per-sample mgatk extraction, using exactly its singlet cells.
process MTDNA_COPY_NUMBER {
  tag "$id"
  cpus 1
  memory '24 GB'
  time '48h'
  publishDir "${params.outdir}/${params.stage}/copy_number", mode: 'copy'
  input:
  tuple val(id), val(arc), val(enriched), val(original), path(metrics), path(json)
  output:
  tuple val(id), path("${id}.mtDNA_copy_number.tsv"), emit: metrics
  tuple val(id), path("${id}.nuclear_reference.bed"), emit: reference
  tuple val(id), path("${id}.settings.txt"), emit: settings
  script:
  def ordinaryArg = original ? "--mgatk-dir " + shellQuote(original) : ''
  def mapFile = "${arc}/per_barcode_metrics.csv"
  def mapArg = file(mapFile).exists() ? "--barcode-map " + shellQuote(mapFile) : ''
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/10_estimate_mtdna_copy_number.R")} --sample ${shellQuote(id)} --bam ${shellQuote(arc + "/atac_possorted_bam.bam")} --fragments ${shellQuote(arc + "/atac_fragments.tsv.gz")} --metrics-tsv ${shellQuote(metrics)} --mgatk-enriched-dir ${shellQuote(enriched)} ${ordinaryArg} ${mapArg} --json ${shellQuote(json)} --out-prefix ${shellQuote(id)}
  """
}

process WNN_STAGE {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path rna
  path atac
  path json
  output:
  path 'wnn.rds'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/12_wnn_stage.R")} wnn.rds ${shellQuote(json)} ${shellQuote(rna)} ${shellQuote(atac)}
  """
}

process ADD_MGATK {
  cpus 2
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path wnn
  path metrics
  path copy_numbers
  output:
  path 'wnn_mgatk.rds'
  path 'all_cells_mtDNA.tsv'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/13_add_mgatk.R")} ${shellQuote(wnn)} wnn_mgatk.rds all_cells_mtDNA.tsv ${metrics.collect { shellQuote(it) }.join(' ')} ${copy_numbers ? "--copy-number-files " + copy_numbers.collect { shellQuote(it) }.join(' ') : ''}
  """
}

process DEG_STAGE {
  cpus 4
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/deg", mode: 'copy'
  input:
  path obj
  output:
  path 'cluster_markers.tsv'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/14_cluster_markers.R")} ${shellQuote(obj)}
  """
}

workflow {
  seenIds = new HashSet()
  // Preserve both mitochondrial sources; step 09 uses enriched, step 10 uses both.
  sample_records = Channel.fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
      if (row.keySet() != ['stage', 'sample_id', 'arc_out', 'mgatk_enriched_dir', 'mgatk_dir'].toSet())
        error 'Sample CSV must contain exactly stage,sample_id,arc_out,mgatk_enriched_dir,mgatk_dir columns'
      row
    }
    .filter { it.stage?.trim() == params.stage }
    .ifEmpty { error "No rows for stage ${params.stage} in ${params.samples}" }
    .map { row ->
      for (key in ['sample_id','arc_out','mgatk_enriched_dir'])
        if (!row[key]?.trim()) error "Missing ${key} in manifest row: ${row}"
      def id = row.sample_id.trim()
      if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error "Invalid sample_id ${id}: start with a letter or number and use only letters, numbers, _, . or -"
      if (!seenIds.add(id)) error "Duplicate sample_id for stage ${params.stage}: ${id}"
      def arc = file(row.arc_out.trim()).toAbsolutePath().toString()
      def enriched = file(row.mgatk_enriched_dir.trim()).toAbsolutePath().toString()
      // Blank ordinary mgatk is allowed: its estimate will be NA until supplied.
      def original = row.mgatk_dir?.trim() ? file(row.mgatk_dir.trim()).toAbsolutePath().toString() : ''
      if (!file("${arc}/filtered_feature_bc_matrix.h5").exists() || !file("${arc}/atac_fragments.tsv.gz").exists())
        error "Missing ARC H5 or fragments: ${arc}"
      if (!file("${arc}/atac_fragments.tsv.gz.tbi").exists()) error "Missing fragment index: ${arc}"
      if (!java.nio.file.Files.isDirectory(file(enriched)))
        error "Missing mgatk_enriched_dir directory: ${enriched}"
      if (params.phase == 'analysis' && !params.copy_number_dir) {
        if (!file("${arc}/atac_possorted_bam.bam").exists())
          error "Step 10 requires ARC atac_possorted_bam.bam: ${arc}"
        if (original && !java.nio.file.Files.isDirectory(file(original)))
          error "Missing mgatk_dir directory: ${original}"
      }
      tuple(id, arc, enriched, original)
    }
    .toSortedList { a, b -> a[0] <=> b[0] }
    .flatMap { it }

  rows = sample_records.map { id, arc, enriched, original -> tuple(id, arc, enriched) }

  PEAK_SAMPLE(rows)
  COMMON_PEAKS(PEAK_SAMPLE.out.map { id, bed -> bed }.collect())
  QC_METRICS(rows.combine(COMMON_PEAKS.out))
  QC_REVIEW(QC_METRICS.out.map { id, tab, plot -> tab }.collect())

  if (params.phase == 'cluster_review') {
    json_review = Channel.value(file(params.json))
    PRE_DOUBLET_REVIEW(rows.combine(COMMON_PEAKS.out).combine(json_review))
  }

  if (params.phase == 'analysis') {
    json = Channel.value(file(params.json))
    reviewDir = params.cluster_review_dir ?: "${params.outdir}/${params.stage}/qc/cluster_review"
    reviewed = rows.combine(COMMON_PEAKS.out).combine(json).map { id, arc, enriched, peaks, cfg ->
      def reviewFile = file("${reviewDir}/${id}.cluster_cells.tsv", checkIfExists: true)
      tuple(id, arc, enriched, peaks, cfg, reviewFile)
    }
    QC_FILTER_DOUBLET(reviewed)
    filtered = QC_FILTER_DOUBLET.out
    MERGE_INTEGRATE_RNA(filtered.map { id, enriched, rna, atac, report, plot -> rna }.collect(), json)
    MERGE_INTEGRATE_ATAC(filtered.map { id, enriched, rna, atac, report, plot -> atac }.collect(), json)
    MGATK_METRICS(filtered.map { id, enriched, rna, atac, report, plot -> tuple(id, enriched, rna) })
    // Use the saved analysis object directly; omitted switches default to true.
    rna_for_wnn = pipelineSettings.RNA_processing?.integrate_rna == false ?
      MERGE_INTEGRATE_RNA.out.unintegrated : MERGE_INTEGRATE_RNA.out.integrated
    atac_for_wnn = pipelineSettings.ATAC_processing?.integrate_atac == false ?
      MERGE_INTEGRATE_ATAC.out.unintegrated : MERGE_INTEGRATE_ATAC.out.integrated
    WNN_STAGE(rna_for_wnn, atac_for_wnn, json)
    // An explicit precomputed directory remains an optional reuse route.
    if (params.copy_number_dir) {
      copy_numbers = rows.map { id, arc, enriched ->
        file("${params.copy_number_dir}/${id}.mtDNA_copy_number.tsv", checkIfExists: true)
      }.collect()
    } else {
      copy_inputs = sample_records.join(MGATK_METRICS.out).combine(json)
      MTDNA_COPY_NUMBER(copy_inputs)
      copy_numbers = MTDNA_COPY_NUMBER.out.metrics.map { id, tab -> tab }.collect()
    }
    ADD_MGATK(WNN_STAGE.out, MGATK_METRICS.out.map { id, tab -> tab }.collect(), copy_numbers)
    DEG_STAGE(ADD_MGATK.out[0])
  }
}
