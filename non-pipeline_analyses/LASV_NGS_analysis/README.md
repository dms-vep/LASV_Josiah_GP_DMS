# Next-generation sequencing analysis pipeline for identifying base calls for Lassa virus sequence

## Contributors

* Caleb Carr

## Description 

A viral deep sequencing analysis pipeline that processes raw sample fastq files and produces a summary of base calls for each sample. This pipeline specifically analyzes the base calling for the Lassa isolate LM395.

## Organization of repository 

* Conda_Envs/ This directory contains all the yaml files that specify the conda environments used for each specific rule.
* Configure/ 
    * Adapter_Files/ This directory contains the default adapter list used to trim adapters automatically.
    * Contaminant_Files/ This directory contains the reference genome for a common next-generation sequencing contaminant (phix).
    * Reference_Files/ This directory contains the text files with the specified urls to download each gene for each reference genome from the NCBI databases.
    * Sample_Files/ This directory contains the csv files that specify the samples to download from the NCBI databases. 
    * __config.yml This yaml file specifies all the pipeline configurations.__ 
* _Data/ This directory is created upon running the pipeline and contains all the reference and sample data._
* Report/ This directory contains all the restructed text files that specify the figure captions in the snakemake report.
* _Results/ This directory is created upon running the pipeline and contains all the processed information and final figures._
* Rules/ This directory contains all the snakemake files that specify the rules for each step of the pipeline.
* Scripts/ This directory contains all the custom python scripts that are used for different steps of the pipeline.
* Software/ This directory contains all the text files and bash scripts needed to create the master conda environment used to run the pipeline.
* Snakefile This snakemake file controls the entire pipeline and specifies the final output of the pipeline.
* cluster.yml This yaml file specifies the cluster configurations used for each rule.
* run_snakemake_cluster.bash This bash file runs the entire pipeline on the cluster.

## Running the pipeline

_Conda is assumed to be installed either through miniconda or anaconda. Please click [here][1] for help to install conda._

[1]: https://docs.anaconda.com/anaconda/install/

### Create master conda environment

```
conda env create -f environment.yml
```

This could fail because the underlying programs may be upgraded creating software conflicts or deprecated commands could result in a non-functioning pipeline. To prevent this, a pinned environment is also included to track the specific program versions.

```
conda env create -f pinned_environment.yml
```

Then activate the conda environment with:

```
conda activate LASV_DeepSeq_Env
```

### Running the pipeline on the Hutch cluster

```shell
sbatch run_snakemake_cluster.bash
```