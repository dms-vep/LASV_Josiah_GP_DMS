# Library design
The code and data used to design the specific mutants in the library would be in [./library_design/](library_design), and are run outside the main pipeline in [Snakefile](Snakefile). However, this study does not contain any specific mutants. 

## Repo organization
 
* `Library_Construction` directory contains directories for library design, initial PCR mutagenesis QC (using Sanger sequencing), and analysis of library titers from transfection and recovery.
* `Experiment_Validations` directory contains scripts, data, plasmid maps, and analysis for validating the effects of mutations to LASV GP, determining antibody IC50s, and test selections. Although some qPCR and flow-based titer data is included for validating selections, etc., the data in this directory typically come from luciferase-based readouts.

## Create DMS library design environment

First, set up the environment used to design the dms library, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [environment.yml](environment.yml).
If you have not previously built the conda environment, then build the environment specified in [environment.yml](environment.yml) with:

    conda env create --file library_design_environment.yml

Then activate it with:

    conda activate lasv_dms_library_design_env

If you've previously built the environment, just do the activation step.