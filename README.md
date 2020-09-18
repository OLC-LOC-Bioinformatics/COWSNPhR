COWSNPhR:
================

Single nucleotide variant detection pipeline using logic modified from https://github.com/USDA-VS/vSNP and SNVPhyl

Documentation is available [here](https://OLC-LOC-Bioinformatics.github.io/COWSNPhR/)

### Quickstart

You must have docker installed, and be a member of the 'docker' users group

Download the DeepVariant image

`docker pull gcr.io/deepvariant-docker/deepvariant:1.0.0`

Install the conda package

`conda install -c olcbioinformatics cowsnphr`

Run the pipeline

`cowsnphr -s /path/to/fastq -r /path/to/ref`

### Reporting issues

If you have any problems installing or running COWSNPhR, or have feature request, please open an issue here on GitHub.
