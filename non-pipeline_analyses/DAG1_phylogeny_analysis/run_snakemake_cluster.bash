#!/bin/bash

# Description:
# bash script to run snakemake pipeline on cluster

# Author:
# Caleb Carr

# Parameters:
# N/A 

# stop on errors
set -e

# Make directory for slurm output
mkdir -p logs_slurm

# Signal that snakemake is running
echo "Running snakemake..."

# Run the main analysis on `slurm` cluster
snakemake \
    -j 99 \
    --use-conda \
    --cluster-config cluster.yml \
    --cluster "sbatch -p {cluster.partition} -c {cluster.cpus} -t {cluster.time} --mem={cluster.mem} -J {cluster.name} -o {cluster.output} -e {cluster.output}" \
    --latency-wait 60 \
    --rerun-triggers mtime 

# Signal that snakemake has complete
echo "Run of snakemake complete."