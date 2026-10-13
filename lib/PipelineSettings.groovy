/** Validate loaded Nextflow parameters, then build only each task's settings. */
class PipelineSettings {
    static final List METRICS = ['nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo',
        'nCount_ATAC', 'nFeature_ATAC', 'nucleosome_signal', 'TSS.enrichment']

    static final List RUN_PARAMETERS = ['samples', 'stage', 'outdir', 'macs3_path',
                                        'conda_env', 'enable_conda', 'phase']
    static final List ANALYSIS_SECTIONS = ['peaks', 'note', 'thresholds', 'sample_overrides',
        'review_QC_cluster', 'doublet', 'RNA_processing', 'ATAC_processing', 'WNN', 'copy_number']

    static Map analysis(Map supplied) {
        object(supplied, 'Pipeline parameters')
        def unknown = supplied.keySet() - (RUN_PARAMETERS + ANALYSIS_SECTIONS)
        if (unknown)
            throw new IllegalArgumentException("Unknown pipeline parameters: ${unknown.join(', ')}")
        // Drop unset sections, preserving nested null QC bounds and copy settings.
        canonical(supplied.subMap(ANALYSIS_SECTIONS).findAll { key, value -> value != null })
    }

    static Map defaults() {
        [review_QC_cluster: [RNA_dims: '1:20', algorithm: 3, resolution: 2],
         doublet: [PC_dims: '1:20', algorithm: 3, resolution: 2, remove_mf_doublets: true],
         RNA_processing: [RNA_dims: '1:20', RNA_umap_dims: '1:50', algorithm: 3,
                          resolution: 2, integrate_rna: true],
         ATAC_processing: [ATAC_dims: '2:40', ATAC_integrate_dims: '1:40', algorithm: 3,
                           resolution: 2, integrate_atac: true, min_cutoff: 'q5'],
         WNN: [RNA_dims: '1:20', ATAC_dims: '2:40', algorithm: 3, resolution: 2,
               plot_point_size: 0.1, lineage_assay: 'SCT',
               lineage_markers: ['Pou3f1', 'Foxa1', 'Pdgfra']],
         copy_number: [bin_size: 500, quantile: 0.985, mapq: 30, baseq: 20, ploidy: 2,
                       min_nuclear_bases: 1000, min_mt_depth: 5, keep_duplicates: false,
                       reference_bed: null, blacklist: null]]
    }

    static void object(value, String label) {
        if (!(value instanceof Map) || value.keySet().any { !(it instanceof String) || !it.trim() })
            throw new IllegalArgumentException("${label} must be an object with nonempty keys")
    }
    static boolean number(value) { value instanceof Number && Double.isFinite(value.doubleValue()) }
    static List dimensions(value, String label) {
        if (!(value instanceof String) || !(value ==~ /[1-9][0-9]*:[1-9][0-9]*/))
            throw new IllegalArgumentException("${label} must be an inclusive range such as 1:20")
        def ends = value.split(':').collect { new BigInteger(it) }
        if (ends[0] > ends[1] || ends[1] > 500)
            throw new IllegalArgumentException("Invalid dimension range: ${label}")
        [ends[0].intValue(), ends[1].intValue()]
    }
    static void bounds(value, String label) {
        object(value, label)
        if (value.keySet() != ['min', 'max'].toSet() || value.values().every { it == null } ||
            !value.values().every { it == null || number(it) } ||
            (value.min != null && value.max != null && value.min > value.max))
            throw new IllegalArgumentException("Invalid min/max bounds: ${label}")
    }
    static void rate(value) {
        if (!(value == 'auto' || (number(value) && value > 0 && value < 1)))
            throw new IllegalArgumentException("Doublet rates must be 'auto' or numbers between zero and one")
    }

    static Map peaks(value) {
        def values = value == null ? [:] : value
        object(values, 'peaks')
        if (values.keySet().any { !(it in ['peak_min', 'peak_max']) })
            throw new IllegalArgumentException('Unknown peaks setting')
        values.each { key, width ->
            if (width != null && !number(width))
                throw new IllegalArgumentException("peaks.${key} must be numeric")
        }
        def selected = [peak_min: values.peak_min == null ? 20 : values.peak_min,
                        peak_max: values.peak_max == null ? 10000 : values.peak_max]
        selected.each { key, width ->
            if (!number(width) || width < 0 || width > Integer.MAX_VALUE ||
                width != Math.floor(width.doubleValue()))
                throw new IllegalArgumentException("${key} must be a nonnegative integer within the R integer range")
            selected[key] = width.intValue()
        }
        if (selected.peak_min >= selected.peak_max)
            throw new IllegalArgumentException('peak_min must be less than peak_max')
        selected
    }

    static Map normalize(Map supplied) {
        supplied = analysis(supplied)
        object(supplied, 'Pipeline JSON')
        object(supplied.thresholds, 'thresholds')
        if (supplied.thresholds.keySet() != METRICS.toSet())
            throw new IllegalArgumentException("thresholds must contain exactly ${METRICS.join(', ')}")
        supplied.thresholds.each { key, value -> bounds(value, "thresholds.${key}") }
        def overrides = supplied.sample_overrides == null ? [:] : supplied.sample_overrides
        object(overrides, 'sample_overrides')
        overrides.each { id, limits ->
            object(limits, "sample_overrides.${id}")
            if (limits.keySet().any { !(it in METRICS) })
                throw new IllegalArgumentException("Unknown QC metric for ${id}")
            limits.each { key, value -> bounds(value, "sample_overrides.${id}.${key}") }
        }
        peaks(supplied.peaks)
        def result = [thresholds: supplied.thresholds, sample_overrides: overrides]
        defaults().each { section, baseline ->
            def settings = supplied[section] == null ? [:] : supplied[section]
            object(settings, section)
            def extra = section == 'doublet' ? ['expected_rate', 'sample_rates', 'exclude_clusters'] : []
            if (settings.keySet().any { !(it in baseline.keySet() || it in extra) })
                throw new IllegalArgumentException("Unknown setting in ${section}")
            result[section] = baseline + settings
        }
        if (result.doublet.PC_dims == null) result.doublet.PC_dims = defaults().doublet.PC_dims
        ['review_QC_cluster', 'doublet', 'RNA_processing', 'ATAC_processing', 'WNN'].each { section ->
            def values = result[section]
            if (!number(values.algorithm) || !(values.algorithm in [1, 2, 3, 4]))
                throw new IllegalArgumentException("${section}.algorithm must be 1, 2, 3, or 4")
            if (!number(values.resolution) || values.resolution <= 0)
                throw new IllegalArgumentException("${section}.resolution must be positive")
            values.each { key, value ->
                if (key.endsWith('_dims')) dimensions(value, "${section}.${key}")
                if (key.startsWith('integrate_') && !(value instanceof Boolean))
                    throw new IllegalArgumentException("${section}.${key} must be true or false")
            }
        }
        def wnnPlot = result.WNN
        if (!number(wnnPlot.plot_point_size) || wnnPlot.plot_point_size <= 0)
            throw new IllegalArgumentException('WNN.plot_point_size must be positive')
        if (!(wnnPlot.lineage_assay instanceof String) || !wnnPlot.lineage_assay.trim())
            throw new IllegalArgumentException('WNN.lineage_assay must be a nonempty assay name')
        def markers = wnnPlot.lineage_markers
        if (!(markers instanceof List) || markers.size() < 1 ||
            !markers.every { it instanceof String && it.trim() } || markers.unique(false).size() != markers.size())
            throw new IllegalArgumentException('WNN.lineage_markers must contain unique nonempty gene symbols')
        def pc = dimensions(result.doublet.PC_dims, 'doublet.PC_dims')
        if (pc[0] != 1 || pc[1] > 50) throw new IllegalArgumentException('doublet.PC_dims must be 1:K with K <= 50')
        def doublet = result.doublet
        doublet.expected_rate = doublet.expected_rate == null ? 'auto' : doublet.expected_rate
        rate(doublet.expected_rate)
        if (!(doublet.remove_mf_doublets instanceof Boolean))
            throw new IllegalArgumentException('doublet.remove_mf_doublets must be true or false')
        ['sample_rates', 'exclude_clusters'].each { key ->
            doublet[key] = doublet[key] == null ? [:] : doublet[key]
            object(doublet[key], "doublet.${key}")
        }
        doublet.sample_rates.values().each { rate(it) }
        doublet.exclude_clusters.values().each { ids ->
            if (!(ids instanceof List) || !ids.every { it instanceof String && it.length() > 0 })
                throw new IllegalArgumentException('exclude_clusters must contain arrays of cluster ID strings')
        }
        def atac = result.ATAC_processing
        def cutoff = atac.min_cutoff
        if (!((number(cutoff) && cutoff >= 0 && cutoff == Math.floor(cutoff.doubleValue())) ||
              (cutoff instanceof String && cutoff ==~ /q[0-9]+([.][0-9]+)?/ &&
               cutoff.substring(1).toBigDecimal() < 100)))
            throw new IllegalArgumentException('ATAC_processing.min_cutoff must be a count or percentile from q0 to below q100')
        if (atac.integrate_atac) {
            def corrected = dimensions(atac.ATAC_integrate_dims, 'ATAC_integrate_dims')
            def used = dimensions(atac.ATAC_dims, 'ATAC_dims')
            def wnn = dimensions(result.WNN.ATAC_dims, 'WNN.ATAC_dims')
            if (corrected[0] != 1 || corrected[1] < Math.max(used[1], wnn[1]))
                throw new IllegalArgumentException('ATAC_integrate_dims must start at 1 and cover ATAC and WNN dimensions')
        }
        result.copy_number.each { key, value ->
            // null means use the standalone script's default, including paths.
            if (value == null) { result.copy_number[key] = defaults().copy_number[key]; return }
            if (key in ['reference_bed', 'blacklist']) {
                if (!(value instanceof String) || !value.startsWith('/'))
                    throw new IllegalArgumentException("copy_number.${key} must be an absolute path or null")
            } else if (key == 'keep_duplicates') {
                if (!(value instanceof Boolean)) throw new IllegalArgumentException('keep_duplicates must be true or false')
            } else {
                if (!number(value) || value < 0 || (key == 'quantile' && value >= 1) ||
                    (key in ['bin_size', 'ploidy'] && value <= 0) ||
                    (key in ['bin_size', 'mapq', 'baseq', 'min_nuclear_bases'] && value != Math.floor(value.doubleValue())))
                    throw new IllegalArgumentException("Invalid copy_number.${key}")
            }
        }
        canonical(result)
    }

    /** Stable key order makes task JSON independent of formatting/key order. */
    static canonical(value) {
        if (value instanceof Map) {
            def result = new LinkedHashMap()
            value.keySet().sort().each { key -> result[key] = canonical(value[key]) }
            return result
        }
        if (value instanceof List) return value.collect { canonical(it) }
        value
    }
    static Map sample(Map settings, String id, boolean doublet = false) {
        def limits = settings.thresholds + (settings.sample_overrides[id] ?: [:])
        if (!doublet) return canonical([process: 'review_QC_cluster', thresholds: limits,
                                       review_QC_cluster: settings.review_QC_cluster])
        def values = settings.doublet
        def selected = values.subMap(['PC_dims', 'algorithm', 'resolution', 'remove_mf_doublets']) +
            [expected_rate: values.sample_rates.containsKey(id) ? values.sample_rates[id] : values.expected_rate,
             exclude_clusters: [(id): values.exclude_clusters[id] ?: []]]
        canonical([process: 'doublet', thresholds: limits, doublet: selected])
    }
    static Map rna(Map settings) {
        def values = settings.RNA_processing
        def count = [dimensions(values.RNA_dims, 'RNA_dims')[1],
                     dimensions(values.RNA_umap_dims, 'RNA_umap_dims')[1],
                     dimensions(settings.WNN.RNA_dims, 'WNN.RNA_dims')[1]].max()
        canonical([process: 'RNA_processing', RNA_processing: values, component_count: count])
    }
    static Map atac(Map settings) {
        def values = new LinkedHashMap(settings.ATAC_processing)
        def count = [dimensions(values.ATAC_dims, 'ATAC_dims')[1],
                     dimensions(settings.WNN.ATAC_dims, 'WNN.ATAC_dims')[1]]
        if (values.integrate_atac) count << dimensions(values.ATAC_integrate_dims, 'ATAC_integrate_dims')[1]
        else values.remove('ATAC_integrate_dims') // unused when integration is disabled
        canonical([process: 'ATAC_processing', ATAC_processing: values, component_count: count.max()])
    }
    static Map wnn(Map settings) {
        canonical([process: 'WNN', WNN: settings.WNN,
                   RNA_processing: [integrate_rna: settings.RNA_processing.integrate_rna],
                   ATAC_processing: [integrate_atac: settings.ATAC_processing.integrate_atac]])
    }
    static Map copyNumber(Map settings) {
        def values = new LinkedHashMap(settings.copy_number)
        if (values.reference_bed) ['bin_size', 'quantile', 'blacklist'].each { values.remove(it) }
        canonical([copy_number: values])
    }
}
