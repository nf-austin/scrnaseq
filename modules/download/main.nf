process DOWNLOAD_REFERENCE {
    tag "Downloading Reference"
    // symlink rather than copy: the 10x reference is multi-GB and duplicating it
    // from work/ into results/ is pure cost on shared storage.
    publishDir "${params.outdir}/reference", mode: 'symlink'

    container 'quay.io/biocontainers/wget:1.25.0'

    input:
    val url

    output:
    path "refdata", emit: ref_dir

    script:
    """
    set -euo pipefail
    mkdir refdata
    wget -q -O - ${url} | tar -xz -C refdata --strip-components=1 || {
        echo "ERROR: could not download the Cell Ranger reference from ${url}" >&2
        echo "       Compute nodes often have no outbound network. Download it once on a" >&2
        echo "       login node and re-run with:  --transcriptome /path/to/refdata" >&2
        exit 1; }
    """

    stub:
    """
    mkdir -p refdata/fasta refdata/genes
    touch refdata/reference.json
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
    wget -q -O probe_set.csv ${url} || {
        echo "ERROR: could not download the probe set from ${url}" >&2
        echo "       Compute nodes often have no outbound network. Download it once on a" >&2
        echo "       login node and re-run with:  --probe_set /path/to/probe_set.csv" >&2
        exit 1; }
    """

    stub:
    """
    echo 'gene_id,probe_seq,probe_id' > probe_set.csv
    """
}
