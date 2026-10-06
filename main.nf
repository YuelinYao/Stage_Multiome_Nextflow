nextflow.enable.dsl = 2

params.samples = null
params.stage = null
params.outdir = 'results'
params.macs3_path = null
params.json = null
params.cluster_review_dir = null
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
if (params.json) {
  if (!file(params.json).exists()) error "Pipeline JSON not found: ${params.json}"
  def settings = new groovy.json.JsonSlurper().parse(file(params.json).toFile())
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
  tuple val(id), val(arc), val(mgatk)
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
  tuple val(id), val(arc), val(mgatk), path(peaks)
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
  path 'pipeline_settings.template.json'
  path 'all_cells_qc_metrics.tsv'
  path 'stage_qc_distributions.pdf'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/04_qc_review.R")} qc_summary.tsv qc_mad_candidates.tsv pipeline_settings.template.json stage_qc_distributions.pdf ${params.peak_min} ${params.peak_max} ${metric_files.collect { shellQuote(it) }.join(' ')}
  """
}

process PRE_DOUBLET_REVIEW {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/cluster_review", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk), path(peaks), path(json)
  output:
  path "${id}.cluster_cells.tsv"
  path "${id}.cluster_summary.tsv"
  path "${id}.cluster_qc.pdf"
  path "${id}.cluster_markers.tsv"
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/05a_review_clusters.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} ${shellQuote(json)}
  """
}

process QC_FILTER_DOUBLET {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/filtered", mode: 'copy', pattern: '*.{qc.tsv,pdf}'
  input:
  tuple val(id), val(arc), val(mgatk), path(peaks), path(json), path(review)
  output:
  tuple val(id), val(mgatk), path("${id}.rna.qc.rds"), path("${id}.atac.qc.rds"), path("${id}.qc.tsv"), path("${id}_Doublet.pdf")
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/05b_qc_doubletfinder.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} ${shellQuote(json)} '${id}.rna.qc.rds' '${id}.atac.qc.rds' '${id}.qc.tsv' ${shellQuote(review)}
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
  path 'rna_stage.rds', emit: stage
  path 'rna_cell_annotations.tsv', emit: annotations
  path 'rna_unintegrated.rds', emit: unintegrated
  path 'rna_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06a_preprocess_rna.R")} ${shellQuote(json)} ${objects.collect { shellQuote(it) }.join(' ')}
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
  path 'atac_stage.rds', emit: stage
  path 'atac_unintegrated.rds', emit: unintegrated
  path 'atac_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06b_preprocess_atac.R")} ${shellQuote(json)} ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MGATK_METRICS {
  tag "$id"
  cpus 2
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/mgatk", mode: 'copy'
  input:
  tuple val(id), val(mgatk), path(qc)
  output:
  path "${id}.mtDNA.tsv"
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06c_mgatk_metrics.R")} ${shellQuote(id)} ${shellQuote(qc)} ${shellQuote(mgatk)} '${id}.mtDNA.tsv'
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
  Rscript ${shellQuote(projectDir.toString() + "/scripts/08_wnn_stage.R")} wnn.rds ${shellQuote(json)} ${shellQuote(rna)} ${shellQuote(atac)}
  """
}

process ADD_MGATK {
  cpus 2
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path wnn
  path metrics
  output:
  path 'wnn_mgatk.rds'
  path 'all_cells_mtDNA.tsv'
  script:
  """
  # R source checksum: ${scriptsHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/09_add_mgatk.R")} ${shellQuote(wnn)} wnn_mgatk.rds all_cells_mtDNA.tsv ${metrics.collect { shellQuote(it) }.join(' ')}
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
  Rscript ${shellQuote(projectDir.toString() + "/scripts/10_cluster_markers.R")} ${shellQuote(obj)}
  """
}

workflow {
  seenIds = new HashSet()
  rows = Channel.fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true)
    .filter { it.stage?.trim() == params.stage }
    .ifEmpty { error "No rows for stage ${params.stage} in ${params.samples}" }
    .map { row ->
      for (key in ['sample_id','arc_out','mgatk_dir'])
        if (!row[key]?.trim()) error "Missing ${key} in manifest row: ${row}"
      def id = row.sample_id.trim()
      if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error "Invalid sample_id ${id}: start with a letter or number and use only letters, numbers, _, . or -"
      if (!seenIds.add(id)) error "Duplicate sample_id for stage ${params.stage}: ${id}"
      def arc = file(row.arc_out.trim()).toAbsolutePath().toString()
      def mt = file(row.mgatk_dir.trim()).toAbsolutePath().toString()
      if (!file("${arc}/filtered_feature_bc_matrix.h5").exists() || !file("${arc}/atac_fragments.tsv.gz").exists())
        error "Missing ARC H5 or fragments: ${arc}"
      if (!file("${arc}/atac_fragments.tsv.gz.tbi").exists()) error "Missing fragment index: ${arc}"
      if (!java.nio.file.Files.isDirectory(file(mt))) error "Missing mgatk directory: ${mt}"
      tuple(id, arc, mt)
    }
    .toSortedList { a, b -> a[0] <=> b[0] }
    .flatMap { it }

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
    reviewed = rows.combine(COMMON_PEAKS.out).combine(json).map { id, arc, mt, peaks, cfg ->
      def reviewFile = file("${reviewDir}/${id}.cluster_cells.tsv", checkIfExists: true)
      tuple(id, arc, mt, peaks, cfg, reviewFile)
    }
    QC_FILTER_DOUBLET(reviewed)
    filtered = QC_FILTER_DOUBLET.out
    MERGE_INTEGRATE_RNA(filtered.map { id, mt, rna, atac, report, plot -> rna }.collect(), json)
    MERGE_INTEGRATE_ATAC(filtered.map { id, mt, rna, atac, report, plot -> atac }.collect(), json)
    MGATK_METRICS(filtered.map { id, mt, rna, atac, report, plot -> tuple(id, mt, rna) })
    WNN_STAGE(MERGE_INTEGRATE_RNA.out.stage, MERGE_INTEGRATE_ATAC.out.stage, json)
    ADD_MGATK(WNN_STAGE.out, MGATK_METRICS.out.collect())
    DEG_STAGE(ADD_MGATK.out[0])
  }
}
