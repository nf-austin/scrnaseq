#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// IMPORTS
// -----------------------------------------------------------------------------
include { DOWNLOAD_REFERENCE; DOWNLOAD_PROBE_SET } from './modules/download'
include { CELLRANGER_COUNT; CELLRANGER_MULTI }     from './modules/cellranger'
include { SCANPY_QC }                              from './modules/scanpy_qc'
include { DUMMY_PROBE_SET }                        from './modules/dummy'

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

// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------
workflow {
    // 1. Ingest inputs
    if (params.run_flex) {
        Channel.fromPath(params.multi_configs)
            | map { file -> tuple(file.baseName, file) }
            | set { ch_samples }
    } else {
        Channel.fromPath(params.fastq_dirs, type: 'dir')
            | map { dir -> tuple(dir.baseName, dir) }
            | set { ch_samples }
    }

    // 2. Resolve references
    if (params.transcriptome) {
        ch_transcriptome = Channel.fromPath(params.transcriptome).first()
    } else {
        ch_transcriptome = DOWNLOAD_REFERENCE(get_transcriptome_url()).first()
    }

    if (params.run_flex) {
        if (params.probe_set) {
            ch_probe_set = Channel.fromPath(params.probe_set).first()
        } else {
            ch_probe_set = DOWNLOAD_PROBE_SET(get_probe_set_url()).first()
        }
    } else {
        ch_probe_set = DUMMY_PROBE_SET().first()
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
                    def i = parts.findIndexOf { it == 'per_sample_outs' }
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