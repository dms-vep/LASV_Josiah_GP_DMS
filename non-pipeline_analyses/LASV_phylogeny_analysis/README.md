# Phylogeny analysis for Lassa virus

This repo analyzes the current Lassa virus sequences that have been sequenced by running a [snakemake](https://snakemake.readthedocs.io/) pipeline that downloads sequences based a list of accessions, processes the sequences, and constructs a phylogenetic tree in addition to other analyses.

Analysis performed by Caleb Carr.

## Genbank accessions

To download the accessions, go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and click *Search by virus*. In the *Search by virus name or taxonomy* box, enter *Mammarenavirus lassaense, taxid:3052310* and hit enter. Then click the  *Download* option, select *Accession List* and *Nucleotide* options and hit *Next*. On the next page, select *Download All Records* and hit *Next*. On the next page, select *Accession with version* and click *Download*. Sequences are downloaded from the list of accessions because more information is extracted from the genbank file during the download process. The current accession list was downloaded on August 10, 2023. 

Additional sequences were downloaded to use as Old World and New World arenavirus outgroups from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and saved in *Outgroup_S_references.fasta* and *Outgroup_S_references_metadata.tsv* files.

```
Outgroup reference S segment genomes:

Mammarenavirus mopeiaense:
NC_006575 from NCBI virus Mopeia viruses ref sequence on 090523
DQ328874 from Amanat et al 2018

Mammarenavirus choriomeningitidis: 
NC_004294 from NCBI virus ref sequence on 090523
NC_077807 from NCBI virus ref sequence on 090523

Mammarenavirus machupoense:
NC_005078 Machupo virus from NCBI virus ref sequence on 090523

Mammarenavirus juninense:
NC_005081 Junin virus from NCBI virus ref sequence on 090523
```

## Snakemake Pipeline

The pipeline can be run automatically using [snakemake](https://snakemake.readthedocs.io/) to run [Snakefile](Snakefile), which reads its configuration from [config.yaml](Configure/config.yml). The results of the automated steps are placed in [Results](Results/).

To run these steps, first build the [conda](https://docs.conda.io/) environment, which installs the necessary programs. First install conda. Then build the environment.

```
conda env create -f environment.yml
```

This could fail because the underlying programs may be upgraded creating software conflicts or deprecated commands could result in a non-functioning pipeline. To prevent this, a pinned environment is also included to track the specific program versions.

```
conda env create -f pinned_environment.yml
```

Then activate the conda environment with:

```
conda activate basic_phylogeny_analysis
```

and then run the pipeline on a computing cluster with [slurm](https://slurm.schedmd.com/documentation.html), which uses the configuration specified in [cluster.yml](cluster.yml):

```
sbatch run_snakemake_cluster.bash
```

## Organization of this repo

- [`Configure`](Configure/): Contains the files needed to configure the pipeline including the [config.yaml](Configure/config.yml) and [Input_Data](Configure/Input_Data/) which contains list of accessions, reference genomes, and outgroup references. 
- [`Rules`](Rules/): Contains the snakemake rules used to run the pipeline.
- [`Scripts`](Scripts/): Contains the custom [python](https://www.python.org/) scripts used for part of the analysis.
- [`Notebooks`](Notebooks/): Contains the custom [Jupyter](https://jupyter.org/) notebooks used for part of the analysis.
- [`Results`](Results/): Contains the automatically generated results from the pipeline analysis.