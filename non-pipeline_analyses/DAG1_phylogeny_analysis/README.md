# Phylogeny analysis of DAG1 gene across all vertebrates

This repo analyzes the current DAG1 vertebrate sequences that have been sequenced by running a [snakemake](https://snakemake.readthedocs.io/) pipeline that downloads sequences based a list of accessions, processes the sequences, and constructs a phylogenetic tree in addition to other analyses.

Analysis performed by Caleb Carr.

## Genbank accessions

To download the accessions, go to [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/). In the *Search* box, enter *DAG1* and hit enter. Then click the  *Orthologs* option. On the *DAG1 - dystroglycan 1* orthologs page, click select all and click *Add to cart* and then click *Download*. In the drop down menu, select *Tabular data (CSV)* and *one sequence per gene* and then click *Download*. Sequences are downloaded from the list of accessions provided in the csv file because more information is extracted from the genbank file during the download process. The current accession list was downloaded on September 29, 2023. 

An additional file was created from just the human and mastomys DAG1 fasta sequences to create a small multiple sequence alignment. Note that the full length DAG1 for mastomys is from *Mastomys coucha* because this is the only fully annotated DAG1 sequence for the Mastomys genus and only alpha-dystroglycan sequences are available for *Mastomys natalensis*. However, it has previously been shown that alpha-dystroglycan is 100% conserved at the amino-acid level across the Mastomys genus [(Tayeh et al)](https://pubmed.ncbi.nlm.nih.gov/20674789/), which was verified by aligning all alpha-dystroglycan protein sequences from Mastomys species. Furthermore, it previously has been shown that only alpha-dystroglycan affects Lassa virus GPC attachment and not beta-dystroglycan [(Kunz et al)](https://pubmed.ncbi.nlm.nih.gov/14644604/).

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
conda activate phylogeny_analysis
```

and then run the pipeline on a computing cluster with [slurm](https://slurm.schedmd.com/documentation.html), which uses the configuration specified in [cluster.yml](cluster.yml):

```
sbatch run_snakemake_cluster.bash
```

## Organization of this repo

- [`Configure`](Configure/): Contains the files needed to configure the pipeline including the [config.yaml](Configure/config.yml) and [Input_Data](Configure/Input_Data/) which contains list of accessions. 
- [`Rules`](Rules/): Contains the snakemake rules used to run the pipeline.
- [`Scripts`](Scripts/): Contains the custom [python](https://www.python.org/) scripts used for part of the analysis.
- [`Notebooks`](Notebooks/): Contains the custom [Jupyter](https://jupyter.org/) notebooks used for part of the analysis.
- [`Results`](Results/): Contains the automatically generated results from the pipeline analysis.