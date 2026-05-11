process SCANPY_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample_id), path(h5_file)
    val run_scrublet
    val species

    output:
    tuple val(sample_id), path("${sample_id}_annotated.h5ad"), emit: h5ad

    script:
    def scrublet_flag = run_scrublet ? "--run_scrublet" : ""
    """
    scanpy_qc.py \\
        --h5 ${h5_file} \\
        --sample_id ${sample_id} \\
        --species ${species} \\
        ${scrublet_flag}
    """
}