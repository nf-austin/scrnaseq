# scanpy stack for SCANPY_QC.
#
# A custom image because no public biocontainer carries scanpy, scrublet,
# leidenalg and igraph together, and SCANPY_QC needs all of them. The Cell Ranger
# and download steps use public images (see nextflow.config).
#
# Before this existed, SCANPY_QC had no container at all and resolved only
# through conda -- which is why the docker/singularity profiles used to set
# conda.enabled = true. Those profiles now disable conda, so this image is what
# makes them correct.
#
# Keep the packages in step with modules/scanpy_qc/environment.yml, the conda
# fallback for the same process.
FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        'scanpy>=1.11' \
        'anndata>=0.11' \
        numpy \
        scipy \
        scikit-learn \
        scrublet \
        h5py \
        leidenalg \
        python-igraph

# Fail the build rather than a run if an import is broken.
RUN python3 -c 'import scanpy, anndata, numpy, scipy, sklearn, h5py, leidenalg; print("All packages OK")'
