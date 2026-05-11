process DUMMY_PROBE_SET {
    output:
    path "dummy_probe_set.csv", emit: probe_set

    script:
    """
    touch dummy_probe_set.csv
    """
}