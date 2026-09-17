#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// IMPORTS
// -----------------------------------------------------------------------------
include { DOWNLOAD_REFERENCE; DOWNLOAD_PROBE_SET } from './modules/download'
include { CELLRANGER_COUNT; CELLRANGER_MULTI }     from './modules/cellranger'
include { SCANPY_QC }                              from './modules/scanpy_qc'

// -----------------------------------------------------------------------------
// FUNCTIONS
// -----------------------------------------------------------------------------
def get_transcriptome_url() {
    if (params.transcriptome_url) return params.transcriptome_url
    if (params.species == 'human') return 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz'
    if (params.species == 'mouse') return 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz'
    error "Unsupported species for auto-download: ${params.species}"
}

def get_probe_set_url() {
    if (params.probe_set_url) return params.probe_set_url

    // Flex v1
    if (params.species == 'human' && params.flex_version == 'v1') return 'https://cf.10xgenomics.com/supp/spatial-exp/probe-set/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'
    if (params.species == 'mouse' && params.flex_version == 'v1') return 'https://cf.10xgenomics.com/supp/spatial-exp/probe-set/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv'

    // Flex v2 (APEX)
    if (params.species == 'human' && params.flex_version == 'v2') return 'https://cf.10xgenomics.com/supp/spatial-exp/probe-set/Chromium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv'
    if (params.species == 'mouse' && params.flex_version == 'v2') return 'https://cf.10xgenomics.com/supp/spatial-exp/probe-set/Chromium_Mouse_Transcriptome_Probe_Set_v2.0_GRCm39-2024-A.csv'

    error "Unsupported combination of species and flex version for auto-download. Please provide --probe_set_url directly."
}


def helpMessage() {
    log.info """
    nf-austin/scrnaseq -- 10x Cell Ranger + species-aware scanpy QC

    Usage, samplesheet (recommended; this is what Seqera Platform launches with):
      nextflow run main.nf -profile docker --input samplesheet.csv --outdir results

      samplesheet.csv columns: sample, fastq_dir, multi_config
      Populate fastq_dir for standard runs, or multi_config with --run_flex true.

    Usage, ad-hoc globs:
      nextflow run main.nf -profile docker --fastq_dirs "data/*"
      nextflow run main.nf -profile docker --run_flex true --multi_configs "configs/*.csv"

    Required (one of):
      --input         Samplesheet CSV.
      --fastq_dirs    Directory glob, one directory of FASTQs per sample.
      --multi_configs CSV glob, one Cell Ranger multi config per run (--run_flex true).

    Common options:
      --outdir        Output directory (default: ${params.outdir}).
      --species       human or mouse (default: ${params.species}).
      --run_flex      Flex / multiplexed mode (default: ${params.run_flex}).
      --run_scrublet  Run Scrublet doublet detection (default: ${params.run_scrublet}).

    On HPC, pre-stage the references -- compute nodes usually have no outbound
    network:
      --transcriptome Pre-downloaded Cell Ranger reference directory.
      --probe_set     Pre-downloaded Flex probe set CSV.

    Add `-profile slurm,singularity --slurm_queue <partition>` and use absolute
    paths. See the README for the HPC notes.
    """.stripIndent()
}

/**
 * Resolve one samplesheet entry to a file or directory.
 *
 * A relative entry is resolved against the samplesheet's OWN directory first,
 * which is what someone editing that sheet expects. Nextflow's default is the
 * launch directory, and on Seqera Platform the launch directory is the work
 * directory -- so a relative path there silently resolves somewhere unrelated.
 * Falls back to launch-dir resolution, and only then reports the entry missing.
 */
def resolveInput(path, sheet_dir, row_num, column) {
    if (path.startsWith('/') || path ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/) {
        return file(path, checkIfExists: true)
    }
    def beside_sheet = sheet_dir.resolve(path)
    if (beside_sheet.exists()) {
        return beside_sheet
    }
    def from_launch = file(path)
    if (from_launch.exists()) {
        return from_launch
    }
    error "Samplesheet row ${row_num}: '${column}' not found as '${beside_sheet}' (relative to the samplesheet) nor as '${from_launch}' (relative to the launch directory). Use an absolute path."
}

/**
 * Turn samplesheet rows into (id, path) tuples.
 *
 * One sheet covers both modes: `fastq_dir` is required when --run_flex is false,
 * `multi_config` when it is true. The unused column may be blank or absent.
 *
 * Validated eagerly over the fully-read row list rather than inside a channel
 * closure: errors raised in a closure are lazy -- they never fire under
 * -preview, and in a real run they surface only once the channel is consumed.
 * A bad samplesheet must fail at launch, before Platform provisions compute.
 */
def buildSamples(rows, sheet_dir) {
    if (!rows) {
        error "Samplesheet is empty: ${params.input}"
    }
    if (!rows[0].containsKey('sample')) {
        error "Samplesheet needs a 'sample' column. Found: ${rows[0].keySet().join(', ')}"
    }
    def col = params.run_flex ? 'multi_config' : 'fastq_dir'
    if (!rows[0].containsKey(col)) {
        error "Samplesheet needs a '${col}' column when --run_flex is ${params.run_flex}. Found: ${rows[0].keySet().join(', ')}"
    }

    def seen = [] as Set
    return rows.withIndex().collect { row, idx ->
        def sample_id = row.sample?.trim()
        if (!sample_id) {
            error "Samplesheet row ${idx + 1} has an empty 'sample' value"
        }
        if (!seen.add(sample_id)) {
            error "Samplesheet has a duplicate sample id: '${sample_id}'. Sample ids become output paths and must be unique."
        }
        def path = row[col]?.trim()
        if (!path) {
            error "Samplesheet row ${idx + 1} ('${sample_id}') has an empty '${col}' value, which is required when --run_flex is ${params.run_flex}"
        }
        tuple(sample_id, resolveInput(path, sheet_dir, idx + 1, col))
    }
}

// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------
workflow {
    if (params.help) {
        helpMessage()
        return
    }

    // 1. Ingest inputs
    def glob_param = params.run_flex ? params.multi_configs : params.fastq_dirs
    def glob_name  = params.run_flex ? '--multi_configs' : '--fastq_dirs'
    if (params.input && glob_param) {
        error "Use either --input (samplesheet) or ${glob_name} (glob), not both."
    }
    if (!params.input && !glob_param) {
        error "No input given. Provide --input samplesheet.csv or ${glob_name}. Run with --help for details."
    }

    if (params.input) {
        def sheet = file(params.input, checkIfExists: true)
        def rows = sheet.splitCsv(header: true, strip: true)
        ch_samples = channel.fromList(buildSamples(rows, sheet.parent))
    }
    else if (params.run_flex) {
        ch_samples = channel.fromPath(params.multi_configs, checkIfExists: true)
            .map { f -> tuple(f.baseName, f) }
    }
    else {
        ch_samples = channel.fromPath(params.fastq_dirs, type: 'dir', checkIfExists: true)
            .map { d -> tuple(d.baseName, d) }
    }

    log.info """
    P I P E L I N E   nf-austin/scrnaseq
    ====================================
    input    : ${params.input ?: glob_param}
    mode     : ${params.run_flex ? 'Flex / multi' : 'standard count'}
    species  : ${params.species}
    scrublet : ${params.run_scrublet}
    outdir   : ${params.outdir}
    """.stripIndent()

    // 2. Resolve references
    if (params.transcriptome) {
        ch_transcriptome = channel.fromPath(params.transcriptome, checkIfExists: true).first()
    } else {
        ch_transcriptome = DOWNLOAD_REFERENCE(get_transcriptome_url()).first()
    }

    // The probe set is only needed in Flex mode. CELLRANGER_COUNT does not take
    // one, so the non-Flex branch defines nothing -- there used to be a
    // DUMMY_PROBE_SET process here whose output was never consumed.
    if (params.run_flex) {
        ch_probe_set = params.probe_set
            ? channel.fromPath(params.probe_set, checkIfExists: true).first()
            : DOWNLOAD_PROBE_SET(get_probe_set_url()).first()
    }

    // 3. Execution routing.
    // CELLRANGER_MULTI can emit multiple per-sample h5s per config, each in a
    // distinct per_sample_outs/<sample>/ directory. Re-key those tuples on the
    // cellranger-internal sample name so downstream QC reports per-sample, not
    // per-config.
    if (params.run_flex) {
        CELLRANGER_MULTI(ch_samples, ch_transcriptome, ch_probe_set)
        CELLRANGER_MULTI.out.h5
            | flatMap { config_id, files ->
                (files instanceof List ? files : [files]).collect { f ->
                    def parts = f.toString().replace('\\', '/').split('/')
                    def i = parts.findIndexOf { part -> part == 'per_sample_outs' }
                    def name = (i >= 0 && i + 1 < parts.size()) ? parts[i + 1] : config_id
                    tuple(name, f)
                }
            }
            | set { ch_h5_flattened }
    } else {
        CELLRANGER_COUNT(ch_samples, ch_transcriptome)
        CELLRANGER_COUNT.out.h5
            | set { ch_h5_flattened }
    }

    // 4. Quality Control
    SCANPY_QC(ch_h5_flattened, params.run_scrublet, params.species)
}