# Submitting Runs to the SRA

Before publishing your study, you need to make the raw sequencing runs avaliable in the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). This directory contains the code and instructions to do this. The recommended way to organize your experiments on the SRA is as follows:

```
BioProject: Variant Library (i.e Omicron BA.1)
|
|___BioSamples: PacBio Barcoding, Individual Studies
	|
	|___SRA Experiments: Individual Runs of Antibodies, Sera, etc..
```

This structure is in keeping with the scheme used for Yeast Display DMS projects like this one [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639956). For example, here is this structure captured in a single run from BioProject to an individual SRA Experiment: 

[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639956) >> [BioSample](https://www.ncbi.nlm.nih.gov/biosample/19925005) >> [SRA Experiment](https://www.ncbi.nlm.nih.gov/sra/SRX11291810[accn])

## Instructions

To upload your samples to the [SRA](https://www.ncbi.nlm.nih.gov/sra), follow along with the code and instructions in the [instructions.ipynb](instructions.ipynb) Jupyter notebook. This notebook will walk you through how to configure this directory, interact with the [SRA](https://www.ncbi.nlm.nih.gov/sra), and prepare your sequencing runs for submission. 