# Deep mutational scanning of Lassa virus Josiah strain glycoprotein using a barcoded lentiviral platform
Study by Kate Crawford, Caleb Carr, and Jesse Bloom.

For documentation of the analysis, see [https://dms-vep.github.io/LASV_Josiah_GP_DMS/](https://dms-vep.github.io/LASV_Josiah_GP_DMS/).

## Organization of this repo

### `dms-vep-pipeline-3` submodule

Most of the analysis is done by the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline-3` submodule.
You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) that you commit back to that repo.

### Code and configuration
The [snakemake](https://snakemake.readthedocs.io/) pipeline itself is run by `dms-vep-pipeline-3/Snakefile` which reads its configuration from [config.yaml](config.yaml).
The [conda](https://docs.conda.io/) environment used by the pipeline is that specified in the `environment.yml` file in [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3).

### Input data
Input data for the pipeline are in [./data/](data).

### Results and documentation
The results of running the pipeline are placed in [./results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in [./docs/](docs), which is rendered via GitHub Pages at [https://dms-vep.github.io/LASV_Josiah_GP_DMS/](https://dms-vep.github.io/LASV_Josiah_GP_DMS/).

### Library design
The design of the mutant library is contained in [./library_design/](library_design).
The design is not part of the pipeline but contains code that must be run separately with its own [conda](https://docs.conda.io/) environment.

## Running the pipeline
To run the pipeline, build the conda environment `dms-vep-pipeline-3` in the `environment.yml` file of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), activate it, and run [snakemake](https://snakemake.readthedocs.io/), such as:

    conda activate dms-vep-pipeline-3
    snakemake -j 32 --use-conda -s dms-vep-pipeline-3/Snakefile

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 32 run_Hutch_cluster.bash