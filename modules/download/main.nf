process DOWNLOAD_REFERENCE {
    tag "Downloading Reference"
    publishDir "${params.outdir}/reference", mode: 'copy'

    container 'quay.io/biocontainers/wget:1.25.0'

    input:
    val url

    output:
    path "refdata", emit: ref_dir

    script:
    """
    set -euo pipefail
    mkdir refdata
    wget -q -O - ${url} | tar -xz -C refdata --strip-components=1
    """
}

process DOWNLOAD_PROBE_SET {
    tag "Downloading Probe Set"
    publishDir "${params.outdir}/reference", mode: 'copy'

    container 'quay.io/nf-core/ubuntu:22.04'

    input:
    val url

    output:
    path "probe_set.csv", emit: probe_set

    script:
    """
    set -euo pipefail
    wget -q -O probe_set.csv ${url}
    """
}