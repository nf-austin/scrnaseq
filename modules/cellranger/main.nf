process CELLRANGER_COUNT {
    tag { sample_id }
    publishDir { "${params.outdir}/cellranger/${sample_id}" }, mode: 'copy'

    container 'quay.io/nf-core/cellranger:10.0.0'

    input:
    tuple val(sample_id), path(fastqs)
    path transcriptome

    output:
    tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix.h5"), emit: h5
    path "${sample_id}/outs/web_summary.html"

    script:
    """
    set -euo pipefail

    # Disable Cell Ranger telemetry (CR >= 9). Ignored on older versions.
    cellranger telemetry disable >/dev/null 2>&1 || true

    # Derive the FASTQ sample prefix from 10x file naming
    # (<prefix>_S<n>_L<lane>_R<read>_001.fastq.gz). Comma-join in the rare
    # case the directory contains multiple prefixes for the same library.
    FASTQ_PREFIX=\$(find -L ${fastqs} -maxdepth 1 -name '*_S*_L*_R*_001.fastq.gz' -printf '%f\\n' \\
        | sed -E 's/_S[0-9]+_L[0-9]+_R[12]_001\\.fastq\\.gz\$//' \\
        | sort -u \\
        | paste -sd, -)

    if [ -z "\${FASTQ_PREFIX}" ]; then
        echo "ERROR: no 10x-style FASTQs found in ${fastqs}" >&2
        exit 1
    fi

    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=${fastqs} \\
        --sample=\${FASTQ_PREFIX} \\
        --expect-cells=${params.expect_cells} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui
    """
}

process CELLRANGER_MULTI {
    tag { sample_id }
    publishDir { "${params.outdir}/cellranger/${sample_id}" }, mode: 'copy'

    container 'nf-core/cellranger:10.0.0'

    input:
    tuple val(sample_id), path(multi_config)
    path transcriptome
    path probe_set

    output:
    tuple val(sample_id), path("${sample_id}/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix.h5"), emit: h5
    path "${sample_id}/outs/per_sample_outs/*/web_summary.html"

    script:
    """
    # Disable Cell Ranger telemetry (CR >= 9). Ignored on older versions.
    cellranger telemetry disable >/dev/null 2>&1 || true

    # Replace placeholders with Nextflow-staged absolute paths
    sed -e "s|REF_PLACEHOLDER|\${PWD}/${transcriptome}|g" \\
        -e "s|PROBE_PLACEHOLDER|\${PWD}/${probe_set}|g" \\
        ${multi_config} > patched_config.csv

    cellranger multi \\
        --id=${sample_id} \\
        --csv=patched_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui
    """
}