process SCANPY_QC {
    tag { sample_id }
    // NOTE: this published path is a contract. nf-austin/echidna discovers
    // scrnaseq output as {scrna_dir}/{sample}/{sample}_annotated.h5ad with
    // --scrna_dir pointing at <outdir>/qc. Do not move or rename it.
    publishDir { "${params.outdir}/qc/${sample_id}" }, mode: 'copy'

    conda "${moduleDir}/environment.yml"
    // Container assigned in nextflow.config's withName block.

    input:
    tuple val(sample_id), path(h5_file)
    val run_scrublet
    val species

    output:
    tuple val(sample_id), path("${sample_id}_annotated.h5ad"), emit: h5ad

    script:
    // Script lives in bin/ and is called bare: Nextflow puts $projectDir/bin on
    // PATH and bind-mounts it into the container.
    def scrublet_flag = run_scrublet ? "--run_scrublet" : ""
    """
    scanpy_qc.py \\
        --h5 ${h5_file} \\
        --sample_id ${sample_id} \\
        --species ${species} \\
        ${scrublet_flag}
    """

    stub:
    """
    touch ${sample_id}_annotated.h5ad
    """
}
