# pub-teddy-human-wgs-telomere

### This repository is develop to estimate the telomere length from the TEDDY whole-genome sequencing (WGS) data using five different tools: Computel, Telseq, Telomerecat, qMotif and Motif_counter.

This repository uses:
- Docker <https://www.docker.com/> to build a container with all the necessary software
- Singularity <https://sylabs.io/docs/> to execute the container in a High Performance Computing environment
- Snakemake <https://snakemake.readthedocs.io/en/stable/> to coordinate the execution of the pipeline
