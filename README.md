# scrnaseq

A Nextflow pipeline for end-to-end single-cell RNA-seq processing: 10x Cell Ranger alignment (standard or Flex/multiplexed chemistries) followed by automated, species-aware QC in scanpy. The pipeline annotates cells with passing QC flags and doublet predictions without dropping raw data, making the resulting `.h5ad` matrices immediately ready for probabilistic modeling in tools like Pyro or scvi-tools.

## Pipeline steps

1. **Reference downloads** — Downloads and decompresses the 10x Genomics reference tarball and/or Flex probe set CSVs based on the `--species` and `--flex_version` (via `wget` inside a small public container). Skipped if local paths are provided.
2. **Alignment & Quantification** (`cellranger count` or `cellranger multi`) — Processes raw FASTQ reads into feature-barcode matrices. For standard runs the FASTQ filename prefix (10x naming `<prefix>_S<n>_L<lane>_R<read>_001.fastq.gz`) is auto-detected from the input directory; the directory basename is only used as the cellranger `--id` / output folder name.
3. **Automated QC** (`scanpy`) — Calculates species-specific mitochondrial/ribosomal/hemoglobin fractions (`MT-`/`RPS`/`RPL`/`HB*` for human, `mt-`/`Rps`/`Rpl`/`Hb*` for mouse), dynamic MAD-based sparsity and overabundance thresholds, and runs doublet detection (`scrublet`). Flags cells via a boolean `passing_qc` mask. Raw counts are preserved.

## Requirements

- Nextflow >= 24.04.0
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
| `--input` | *(recommended)* | Samplesheet CSV with `sample`, `fastq_dir`, `multi_config`. |
| `--slurm_queue` | *(cluster default)* | SLURM partition. Used by `-profile slurm`. |
| `--slurm_account` | *(none)* | SLURM account to charge. |
| `--cluster_options` | *(none)* | Raw sbatch options; use the `=` form for `--`-prefixed values. |
| `--singularity_cache_dir` | `$NXF_SINGULARITY_CACHEDIR` | Shared directory for pulled images. |
| `--conda_cache_dir` | `$NXF_CONDA_CACHEDIR` | Shared directory for conda environments. |
| `--singularity_bind` | *(none)* | Extra bind mounts, e.g. `/mnt/gpfs`. |
| `--scanpy_container` | `ghcr.io/nf-austin/scrnaseq:1.0.0` | Image for `SCANPY_QC`. |


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
| `--create_bam` | `true` | Whether `cellranger count` writes BAM output. Required by Cell Ranger >= 8 (no default in the tool itself). Set `false` to skip BAM generation for faster, smaller runs. |

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

## Seqera Platform (Nextflow Tower)

The repo ships everything Platform needs:

- **`nextflow_schema.json`** — renders the launch form. `--input` appears as a file picker wired to
  Data Explorer, options are grouped by stage, and tuning knobs are marked hidden.
- **`assets/schema_input.json`** — the samplesheet contract (`sample`, `fastq_dir`, `multi_config`).
  Two example sheets ship in `assets/`: one standard, one Flex.
- **`tower.yml`** — puts the Cell Ranger web summaries, the QC-annotated h5ads and the Nextflow
  execution report in the run's **Reports** tab.

To add it: **Pipelines → Add pipeline**, point at this repository, and pick a compute environment.
Use **absolute paths** for `--input`, the FASTQ directories it references, `--transcriptome` and
`--outdir`.

### A note on `sample` in Flex mode

In standard mode the sheet's `sample` is the sample id, and it is what appears in `results/`. **In
Flex mode it names the *config*.** The per-sample outputs are keyed by the sample names Cell Ranger
itself reports under `per_sample_outs/`, which come from inside the multi config — so the names in
`results/qc/` are those, not necessarily the `sample` column. This is a Cell Ranger behaviour the
samplesheet cannot override.

## HPC / SLURM

The `slurm` profile sets only the executor and queue, so it composes with an engine profile in
either order:

```bash
nextflow run nf-austin/scrnaseq \
    -profile slurm,singularity \
    --slurm_queue normal \
    --input /mnt/gpfs/project/sheet.csv \
    --transcriptome /mnt/gpfs/refs/refdata-gex-GRCh38-2024-A \
    --outdir /mnt/gpfs/project/results \
    --singularity_cache_dir /mnt/gpfs/shared/singularity
```

- **Pass `--transcriptome` (and `--probe_set` for Flex).** `DOWNLOAD_REFERENCE` fetches a multi-GB
  10x reference *inside a task*, and compute nodes on most clusters have no outbound network. Fetch
  it once on a login node; both download steps fail with a message naming the flag to use.
- **Put `--singularity_cache_dir` on shared storage.** The Cell Ranger image alone is multi-GB and
  `$HOME` is usually quota-limited and not always mounted on compute nodes.
- **`--singularity_bind` matters more here than elsewhere.** `CELLRANGER_MULTI` writes absolute host
  paths into its patched config via `sed`, so on a symlinked filesystem (`/data` → `/mnt/gpfs/...`)
  the container resolves them to nothing. Bind the real parent: `--singularity_bind /mnt/gpfs`.
- **Quote option values that start with `--` using the `=` form**, e.g.
  `--cluster_options='--qos=long'`. The space form is parsed by Nextflow as a bare flag.
- **Seqera Platform already sets the executor** when you launch against a SLURM compute environment,
  so `-profile slurm` is mainly for launching by hand from a login node.

## Container images

| Image | Source | Used by |
| --- | --- | --- |
| `quay.io/nf-core/cellranger:10.0.0` | public | `CELLRANGER_COUNT`, `CELLRANGER_MULTI` |
| `quay.io/biocontainers/wget:1.25.0` | public | `DOWNLOAD_REFERENCE` |
| `quay.io/nf-core/ubuntu:22.04` | public | `DOWNLOAD_PROBE_SET` |
| `ghcr.io/nf-austin/scrnaseq:<ver>` | `modules/scanpy_qc/Dockerfile` | `SCANPY_QC` |

Only `SCANPY_QC` needs a built image — no public biocontainer carries scanpy, scrublet, leidenalg
and igraph together. Everything else uses a public image. The GHCR package must be **public** for
`nextflow run` to pull it without credentials.

**`-profile conda` cannot run this pipeline end to end**: Cell Ranger is proprietary and is only
distributed as a container image. `SCANPY_QC` does ship an `environment.yml`.

## Notes

- `nextflow run . -stub-run --input assets/samplesheet_example.csv --transcriptome <dir>` exercises
  the real channel wiring and publishing with no containers and no data.
- `nextflow lint main.nf nextflow.config modules/*/main.nf` catches config errors that `-preview`
  accepts.
- **Downstream contract:** `nf-austin/echidna` discovers this pipeline's output as
  `{scrna_dir}/{sample}/{sample}_annotated.h5ad` with `--scrna_dir <outdir>/qc`. Do not move or
  rename `SCANPY_QC`'s published path.
