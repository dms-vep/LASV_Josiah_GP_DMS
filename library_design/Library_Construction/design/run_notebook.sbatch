#!/bin/bash
#SBATCH
#SBATCH -o lib_design.out
#SBATCH -e lib_design.err
#SBATCH -c 4
#SBATCH -t 24:00:00

runipy -o design.ipynb
