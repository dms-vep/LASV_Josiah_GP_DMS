# Library construction

This directory contains sub-directories for the steps of making the LASV GP
pseudotyped lentivirus DMS library.

## Directories:

### `design`:
    
    Contains:
    1. Notebooks, alignments, and trees for initial analysis of Arenavirus GP diversity used to determine what regions of LASV GP to mutate.
    Decided to not exclude any parts of GP.

    2. Directory containing the script, sequence, and primers for the final library design.

### `sanger_analysis`:
    Contains:
    1. Lists of mutations for each round of mutagenesis as well as the script and resulting analyses of the results of Sanger sequencing >8 clones per library.
        Kate's local version contains the `.ab1` trace files, but these are not tracked by GitHub.

## Create DMS library design environment

First, set up the environment used to design the dms library, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [library_design_environment.yaml](library_design_environment.yaml).
If you have not previously built the conda environment, then build the environment specified in [library_design_environment.yaml](library_design_environment.yaml) with:

    conda env create --file library_design_environment.yml

Then activate it with:

    conda activate lasv_dms_library_design_env

If you've previously built the environment, just do the activation step.



