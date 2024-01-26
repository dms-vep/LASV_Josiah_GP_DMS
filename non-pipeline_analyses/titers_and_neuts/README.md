# Pseudovirus titers and neutralization assays

This repo analyzes a combination of flow-based and luciferase-based results and creates figures.

Analysis performed by Caleb Carr.

## Activate conda environment

To run the notebook, first build the [conda](https://docs.conda.io/) environment, which installs the necessary programs. First install conda. Then build the environment.

```
conda env create -f virus_titers_and_neuts.yml
```

Then activate the conda environment with:

```
conda activate virus_titers_and_neuts
```

## Organization of this repo

- [`data`](data/): Contains the raw pseudovirus titers and fraction infectivity values.
- [`results`](results/): Contains the output figures of titers and neutralization assays.