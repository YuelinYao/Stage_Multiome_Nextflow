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
// Serialize only this task's val settings. The supplied pipeline JSON is never staged.
def taskSettingsFile(settings) {
  "cat > process_settings.json <<'PROCESS_SETTINGS_JSON'\n" +
    groovy.json.JsonOutput.toJson(settings) + "\nPROCESS_SETTINGS_JSON"
}
// Track shared files as path inputs, while R keeps their original absolute paths.
// Staging creates symlinks; BAMs/fragments are not copied into task directories.
def arcDataFiles(arc, boolean summary = false) {
  def names = ['filtered_feature_bc_matrix.h5', 'atac_fragments.tsv.gz', 'atac_fragments.tsv.gz.tbi']
  if (summary) names << 'summary.csv'
  names.collect { file("${arc}/${it}", checkIfExists: true) }
}
def mgatkDataFiles(directory, boolean depthOnly = false) {
  def names = file(directory).toFile().listFiles()?.findAll { it.isFile() } ?: []
  def patterns = depthOnly ? [/[.]depthTable[.]txt$/] :
    [/[.]A[.]txt[.]gz$/, /[.]C[.]txt[.]gz$/, /[.]G[.]txt[.]gz$/, /[.]T[.]txt[.]gz$/,
     /_refAllele[.]txt(?:[.]gz)?$/, /[.]depthTable[.]txt$/]
  patterns.collect { pattern ->
    def matches = names.findAll { it.name =~ pattern }
    if (matches.size() != 1)
      error "Expected exactly one mgatk file matching ${pattern} in ${directory}"
    file(matches[0].toPath())
  }.sort { it.name }
}
def copyNumberDataFiles(arc, enriched, original, settings) {
  def paths = [file("${arc}/atac_possorted_bam.bam", checkIfExists: true)]
  paths.addAll(mgatkDataFiles(enriched, true))
  if (original) paths.addAll(mgatkDataFiles(original, true))
  def map = file("${arc}/per_barcode_metrics.csv")
  if (map.exists()) paths << map
  def copy = settings.copy_number
  if (copy.reference_bed) paths << file(copy.reference_bed, checkIfExists: true)
  else {
    paths << file("${arc}/atac_fragments.tsv.gz", checkIfExists: true)
    if (copy.blacklist) paths << file(copy.blacklist, checkIfExists: true)
  }
  paths.unique().sort { it.toString() }
}
// Include R sources in task scripts so -resume notices helper/script edits.
def scriptsHash() {
  def sourceDigest = java.security.MessageDigest.getInstance('SHA-256')
  file("${projectDir}/scripts").toFile().listFiles()
    .findAll { it.name.endsWith('.R') }.sort { it.name }.each {
      sourceDigest.update(it.name.getBytes('UTF-8'))
      sourceDigest.update(it.bytes)
    }
  return sourceDigest.digest().encodeHex().toString()

}

process PEAK_SAMPLE {
  tag "$id"
  cpus 4
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/peaks/sample", mode: 'copy'
  input:
  tuple val(id), val(arc), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.bed")
  script:
  """
  # R source checksum: ${scriptsHash()}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/01_peaks.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(params.macs3_path)} '${id}.bed'
  """
}

process COMMON_PEAKS {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/peaks", mode: 'copy'
  input:
  path beds
  val peakMin
  val peakMax
  output:
  path 'common_peaks.bed'
  script:
  """
  # R source checksum: ${scriptsHash()}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/02_consensus_peaks.R")} common_peaks.bed ${peakMin} ${peakMax} ${beds.collect { shellQuote(it) }.join(' ')}
  """
}

process QC_METRICS {
  tag "$id"
  cpus 4
  memory '48 GB'
  input:
  tuple val(id), val(arc), path(peaks), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.metrics.tsv")
  script:
  """
  # R source checksum: ${scriptsHash()}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/03_qc_metrics.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} '${id}.metrics.tsv'
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
  # R source checksum: ${scriptsHash()}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/04_qc_review.R")} qc_summary.tsv qc_mad_candidates.tsv stage_qc_distributions.pdf ${metric_files.collect { shellQuote(it) }.join(' ')}
  """
}

process PRE_DOUBLET_REVIEW {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/cluster_review", mode: 'copy'
  input:
  tuple val(id), val(arc), path(peaks), val(settings), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  path "${id}.cluster_cells.tsv"
  path "${id}.cluster_summary.tsv"
  path "${id}.cluster_qc.pdf"
  path "${id}.cluster_markers.tsv"
  script:
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/05_review_clusters.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} process_settings.json
  """
}

process QC_FILTER_DOUBLET {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/filtered", mode: 'copy', pattern: '*.{qc.tsv,pdf}'
  input:
  tuple val(id), val(arc), path(peaks), val(settings), path(review), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.rna.qc.rds"), path("${id}.atac.qc.rds"), path("${id}.qc.tsv"), path("${id}_Doublet.pdf")
  script:
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06_qc_doubletfinder.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} process_settings.json '${id}.rna.qc.rds' '${id}.atac.qc.rds' '${id}.qc.tsv' ${shellQuote(review)}
  """
}

process MERGE_INTEGRATE_RNA {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'rna_*'
  input:
  path objects
  val settings
  output:
  path 'rna_cell_annotations.tsv', emit: annotations
  path 'rna_unintegrated.rds', emit: unintegrated
  path 'rna_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/07_preprocess_rna.R")} process_settings.json ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MERGE_INTEGRATE_ATAC {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'atac_*.rds'
  input:
  path objects
  val settings
  output:
  path 'atac_unintegrated.rds', emit: unintegrated
  path 'atac_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/08_preprocess_atac.R")} process_settings.json ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MGATK_METRICS {
  tag "$id"
  cpus 2
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/mgatk", mode: 'copy'
  input:
  tuple val(id), val(mgatk_enriched), path(qc), path(mgatk_inputs, stageAs: 'mgatk_inputs/*')
  output:
  tuple val(id), path("${id}.mtDNA.tsv")
  script:
  """
  # R source checksum: ${scriptsHash()}
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
  tuple val(id), val(arc), val(enriched), val(original), path(metrics), val(settings), path(data_inputs, stageAs: 'external_inputs/input??/*')
  output:
  tuple val(id), path("${id}.mtDNA_copy_number.tsv"), emit: metrics
  tuple val(id), path("${id}.nuclear_reference.bed"), emit: reference
  tuple val(id), path("${id}.settings.txt"), emit: settings
  script:
  def ordinaryArg = original ? "--mgatk-dir " + shellQuote(original) : ''
  def mapFile = "${arc}/per_barcode_metrics.csv"
  def mapArg = file(mapFile).exists() ? "--barcode-map " + shellQuote(mapFile) : ''
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/10_estimate_mtdna_copy_number.R")} --sample ${shellQuote(id)} --bam ${shellQuote(arc + "/atac_possorted_bam.bam")} --fragments ${shellQuote(arc + "/atac_fragments.tsv.gz")} --metrics-tsv ${shellQuote(metrics)} --mgatk-enriched-dir ${shellQuote(enriched)} ${ordinaryArg} ${mapArg} --json process_settings.json --out-prefix ${shellQuote(id)}
  """
}

process WNN_STAGE {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path rna
  path atac
  val settings
  output:
  path 'wnn.rds'
  script:
  """
  # R source checksum: ${scriptsHash()}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/12_wnn_stage.R")} wnn.rds process_settings.json ${shellQuote(rna)} ${shellQuote(atac)}
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
  # R source checksum: ${scriptsHash()}
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
  # R source checksum: ${scriptsHash()}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/14_cluster_markers.R")} ${shellQuote(obj)}
  """
}

workflow {
  if (!params.samples || !params.stage) error 'Provide --samples manifest.csv and --stage E8.5'
  if (!(params.phase in ['review', 'cluster_review', 'analysis'])) error '--phase must be review, cluster_review or analysis'
  if (!(params.stage.toString() ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Stage must be a safe directory name'
  // Read only peak settings here: a review JSON can still have unset QC bounds.
  def peakSettings = [:]
  def pipelineSettings = [:]
  if (params.json) {
    if (!file(params.json).exists()) error "Pipeline JSON not found: ${params.json}"
    def settings = PipelineSettings.read(file(params.json).toFile())
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
  def selectedPeaks = PipelineSettings.peaks(peakSettings, params.peak_min, params.peak_max)
  def peakMin = selectedPeaks.peak_min
  def peakMax = selectedPeaks.peak_max
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

  // Validate the complete supplied config once before projecting task inputs.
  // Initial review intentionally permits unfinished QC bounds.
  def scopedSettings = params.phase == 'review' ? [:] : PipelineSettings.normalize(pipelineSettings)
  def rnaSettings = params.phase == 'analysis' ? PipelineSettings.rna(scopedSettings) : [:]
  def atacSettings = params.phase == 'analysis' ? PipelineSettings.atac(scopedSettings) : [:]
  def wnnSettings = params.phase == 'analysis' ? PipelineSettings.wnn(scopedSettings) : [:]
  def copyNumberSettings = params.phase == 'analysis' ? PipelineSettings.copyNumber(scopedSettings) : [:]

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
      ['sample_id', 'arc_out', 'mgatk_enriched_dir'].each { key ->
        if (!row[key]?.trim()) error "Missing ${key} in manifest row: ${row}"
      }
      def id = row.sample_id.trim()
      if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error "Invalid sample_id ${id}: start with a letter or number and use only letters, numbers, _, . or -"
      if (id == 'ALL_STAGE') error 'ALL_STAGE is reserved for the combined QC summary'
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
      if (params.phase == 'analysis') {
        if (!file("${arc}/summary.csv").exists()) error "Missing ARC summary.csv: ${arc}"
        mgatkDataFiles(enriched)
      }
      if (params.phase == 'analysis' && !params.copy_number_dir) {
        if (!file("${arc}/atac_possorted_bam.bam").exists())
          error "Step 10 requires ARC atac_possorted_bam.bam: ${arc}"
        if (original && !java.nio.file.Files.isDirectory(file(original)))
          error "Missing mgatk_dir directory: ${original}"
        copyNumberDataFiles(arc, enriched, original, copyNumberSettings)
      }
      tuple(id, arc, enriched, original)
    }
    .toSortedList { a, b -> a[0] <=> b[0] }
    .flatMap { it }

  rows = sample_records.map { id, arc, enriched, original -> tuple(id, arc, enriched) }

  peak_inputs = rows.map { id, arc, enriched -> tuple(id, arc, arcDataFiles(arc)) }
  PEAK_SAMPLE(peak_inputs)
  COMMON_PEAKS(PEAK_SAMPLE.out.map { id, bed -> bed }.collect().map { it.sort { path -> path.name } }, peakMin, peakMax)
  metric_inputs = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
    tuple(id, arc, peaks, arcDataFiles(arc))
  }
  QC_METRICS(metric_inputs)
  QC_REVIEW(QC_METRICS.out.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } })

  if (params.phase == 'cluster_review') {
    review_inputs = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
      tuple(id, arc, peaks, PipelineSettings.sample(scopedSettings, id), arcDataFiles(arc))
    }
    PRE_DOUBLET_REVIEW(review_inputs)
  }

  if (params.phase == 'analysis') {
    reviewDir = params.cluster_review_dir ?: "${params.outdir}/${params.stage}/qc/cluster_review"
    reviewed = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
      def reviewFile = file("${reviewDir}/${id}.cluster_cells.tsv", checkIfExists: true)
      tuple(id, arc, peaks, PipelineSettings.sample(scopedSettings, id, true), reviewFile, arcDataFiles(arc, true))
    }
    QC_FILTER_DOUBLET(reviewed)
    filtered = QC_FILTER_DOUBLET.out
    MERGE_INTEGRATE_RNA(filtered.map { id, rna, atac, report, plot -> rna }.collect().map { it.sort { path -> path.name } }, rnaSettings)
    MERGE_INTEGRATE_ATAC(filtered.map { id, rna, atac, report, plot -> atac }.collect().map { it.sort { path -> path.name } }, atacSettings)
    mgatk_inputs = sample_records.join(filtered).map { id, arc, enriched, original, rna, atac, report, plot ->
      tuple(id, enriched, rna, mgatkDataFiles(enriched))
    }
    MGATK_METRICS(mgatk_inputs)
    // Use the saved analysis object directly; omitted switches default to true.
    rna_for_wnn = scopedSettings.RNA_processing.integrate_rna == false ?
      MERGE_INTEGRATE_RNA.out.unintegrated : MERGE_INTEGRATE_RNA.out.integrated
    atac_for_wnn = scopedSettings.ATAC_processing.integrate_atac == false ?
      MERGE_INTEGRATE_ATAC.out.unintegrated : MERGE_INTEGRATE_ATAC.out.integrated
    WNN_STAGE(rna_for_wnn, atac_for_wnn, wnnSettings)
    // An explicit precomputed directory remains an optional reuse route.
    if (params.copy_number_dir) {
      copy_numbers = rows.map { id, arc, enriched ->
        file("${params.copy_number_dir}/${id}.mtDNA_copy_number.tsv", checkIfExists: true)
      }.collect().map { it.sort { path -> path.name } }
    } else {
      copy_inputs = sample_records.join(MGATK_METRICS.out).map { id, arc, enriched, original, metrics ->
        tuple(id, arc, enriched, original, metrics, copyNumberSettings,
          copyNumberDataFiles(arc, enriched, original, copyNumberSettings))
      }
      MTDNA_COPY_NUMBER(copy_inputs)
      copy_numbers = MTDNA_COPY_NUMBER.out.metrics.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } }
    }
    ADD_MGATK(WNN_STAGE.out, MGATK_METRICS.out.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } }, copy_numbers)
    DEG_STAGE(ADD_MGATK.out[0])
  }
}
