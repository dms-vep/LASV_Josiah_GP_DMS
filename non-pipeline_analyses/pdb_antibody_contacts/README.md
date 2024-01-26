# Antibody distance to GPC analysis

This repo analyzes the distance between antibodies and GPC.

Analysis performed by Caleb Carr.

## Activate conda environment

To run the notebook, first build the [conda](https://docs.conda.io/) environment, which installs the necessary programs. First install conda. Then build the environment.

```
conda env create -f environment.yml
```

Then activate the conda environment with:

```
conda activate pdb_r_analysis
```

## Organization of this repo

- [`data`](data/): Contains the PDB files for antibodies and GPC structures.
- [`results`](results/): Contains the output files with antibody contact sites and distances.