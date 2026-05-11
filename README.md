# scrnaseq

A Nextflow pipeline for end-to-end single-cell RNA-seq processing: 10x Cell Ranger alignment (standard or Flex/multiplexed chemistries) followed by automated, species-aware QC in scanpy. The pipeline annotates cells with passing QC flags and doublet predictions without dropping raw data, making the resulting `.h5ad` matrices immediately ready for probabilistic modeling in tools like Pyro or scvi-tools.

## Pipeline steps

1. **Reference downloads** — Downloads and decompresses the 10x Genomics reference tarball and/or Flex probe set CSVs based on the `--species` and `--flex_version` (via `curl` inside an Ubuntu container). Skipped if local paths are provided.
2. **Alignment & Quantification** (`cellranger count` or `cellranger multi`) — Processes raw FASTQ reads into feature-barcode matrices. For standard runs the FASTQ filename prefix (10x naming `<prefix>_S<n>_L<lane>_R<read>_001.fastq.gz`) is auto-detected from the input directory; the directory basename is only used as the cellranger `--id` / output folder name.
3. **Automated QC** (`scanpy`) — Calculates species-specific mitochondrial/ribosomal/hemoglobin fractions (`MT-`/`RPS`/`RPL`/`HB*` for human, `mt-`/`Rps`/`Rpl`/`Hb*` for mouse), dynamic MAD-based sparsity and overabundance thresholds, and runs doublet detection (`scrublet`). Flags cells via a boolean `passing_qc` mask. Raw counts are preserved.

## Requirements

- Nextflow >= 23.04
- Docker or Singularity (Cell Ranger runs in the `nf-core/cellranger:10.0.0` container)
- Conda / Mamba (the scanpy QC step is provisioned from `modules/scanpy_qc/environment.yml`; conda is auto-enabled by the `docker` and `singularity` profiles, or use `-profile conda` standalone)

Cell Ranger telemetry is disabled at the start of every `cellranger` invocation (`cellranger telemetry disable`), with a no-op fallback for versions that predate the subcommand.

## Usage

You can run the pipeline either from a local checkout (`main.nf`) or directly from GitHub once a release tag is published:

```bash
nextflow run nf-austin/scrnaseq -r v1.0.0 ...
```

### Standard scRNA-seq (Human, auto-download reference)

```bash
nextflow run main.nf \
    -profile docker \
    --fastq_dirs "data/*" \
    --species human
```

Notes:
- Each `data/*` directory should contain FASTQs for a single library, named per the 10x convention (`<prefix>_S<n>_L<lane>_R<read>_001.fastq.gz`). The FASTQ prefix is auto-detected.
- The directory basename becomes the cellranger `--id` and the QC `sample_id` used downstream.

### Mouse data

```bash
nextflow run main.nf \
    -profile docker \
    --fastq_dirs "data/*" \
    --species mouse
```

### Flex / Multiplexed mode (v1 or v2/APEX)

To run Flex data, you must provide `.csv` configuration files for `cellranger multi`.

Inside your CSV files, use `REF_PLACEHOLDER` for the reference path and `PROBE_PLACEHOLDER` for the probe-set path. The pipeline automatically overwrites these with the correct Nextflow-staged paths during execution.

Example `configs/sample_1.csv`:

```
[gene-expression]
reference,REF_PLACEHOLDER
probe-set,PROBE_PLACEHOLDER

[libraries]
fastq_id,fastqs,feature_types
sample_1,/path/to/fastqs,Multiplexing Capture
```

Run APEX (Flex v2) for Human:

```bash
nextflow run main.nf \
    -profile docker \
    --run_flex true \
    --multi_configs "configs/*.csv" \
    --species human \
    --flex_version v2
```

In Flex mode each config's basename is the cellranger `--id` (the output folder under `results/cellranger/`), but a single `cellranger multi` run can emit multiple per-sample matrices (one per `per_sample_outs/<sample>/`). QC is run once per per-sample matrix, and each is keyed on the cellranger-internal sample name rather than the config basename.

### Supply pre-downloaded references

```bash
nextflow run main.nf \
    -profile docker \
    --fastq_dirs "data/*" \
    --transcriptome "/path/to/refdata-gex-GRCh38-2024-A" \
    --probe_set "/path/to/probe_set.csv"
```

### Toggle Scrublet

Scrublet is run by default. To bypass doublet detection:

```bash
nextflow run main.nf \
    -profile docker \
    --fastq_dirs "data/*" \
    --run_scrublet false
```

## Parameters

### Core inputs

| Parameter | Default | Description |
| --- | --- | --- |
| `--fastq_dirs` | `data/*` | Glob pattern for sample directories containing FASTQs (standard mode). |
| `--multi_configs` | `configs/*.csv` | Glob pattern for Cell Ranger multi configuration files (flex mode). |
| `--outdir` | `results` | Output directory. |

### Biology & chemistry

| Parameter | Default | Description |
| --- | --- | --- |
| `--species` | `human` | Species target. Options: `human`, `mouse`. Updates QC gene parsing and reference defaults. |
| `--flex_version` | `v1` | Flex chemistry version. Options: `v1`, `v2` (APEX). Updates probe set URL defaults. |

### References (auto-resolved if left null)

| Parameter | Default | Description |
| --- | --- | --- |
| `--transcriptome` | `null` | Path to a local Cell Ranger reference directory. |
| `--probe_set` | `null` | Path to a local Flex probe set CSV file. |
| `--transcriptome_url` | `null` | Override the default download URL for the transcriptome. |
| `--probe_set_url` | `null` | Override the default download URL for the probe set CSV. |

### Pipeline mode

| Parameter | Default | Description |
| --- | --- | --- |
| `--run_flex` | `false` | Enable `cellranger multi` mode. |
| `--run_scrublet` | `true` | Execute Scrublet for doublet prediction during the QC process. |
| `--expect_cells` | `10000` | Target cell recovery for `cellranger count`. |

### Resource limits

| Parameter | Default | Description |
| --- | --- | --- |
| `--max_memory` | `128.GB` | Maximum memory available to any single process. |
| `--max_cpus` | `32` | Maximum CPUs available to any single process. |
| `--max_time` | `72.h` | Maximum wall time available to any single process. |

## Output structure

```
results/
├── reference/                # downloaded references and probe sets
├── cellranger/               # per-id cellranger outputs (one folder per --id)
│   └── {id}/
│       └── outs/             # standard: filtered_feature_bc_matrix.h5, web_summary.html
│                             # flex:     per_sample_outs/<sample>/{count, web_summary.html}
└── qc/                       # per-sample scanpy qc outputs
    └── {sample_id}/
        └── {sample_id}_annotated.h5ad   # AnnData with passing_qc + doublet flags
```

For standard runs `{id}` and `{sample_id}` are the same (the FASTQ directory basename). For Flex runs `{id}` is the multi-config basename and `{sample_id}` is the per-sample name emitted by `cellranger multi`.