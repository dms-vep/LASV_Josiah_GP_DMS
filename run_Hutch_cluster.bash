#!/bin/bash
#
#SBATCH -c 32

snakemake -j 32 --use-conda --conda-frontend conda --rerun-incomplete -s dms-vep-pipeline-3/Snakefile