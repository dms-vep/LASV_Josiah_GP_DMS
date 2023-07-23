# Input data

## PacBio full-length variant sequencing
[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file having features to parse with [alignparse](https://jbloomlab.github.io/alignparse/). Must have *gene* (the gene of interest) and *barcode* features.

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/).

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants.
It must have the following columns:
 - `library`: name of the library sequenced
    + *LibA*: concatenated A1-48, A2-1, and A2-2 because these represented three sorts of the same pool of cells
    + *LibB*: concatenated B1-48, B2-1, and B2-2 because these represented three sorts of the same pool of cells
 - `run`: date of the pacbio library submission (use this date to refer to experimental notebook)
    + *A1-48*: 210423
    + *A2-1*: 210430
    + *A2-2*: 210430
    + *B1-48*: 210423
    + *B2-1*: 210430
    + *B2-2*: 210430
    + *LibA*: concatenated on 220404
    + *LibB*: concatenated on 220404
 - `fastq`: FASTQ file from running CCS
    + Original PacBio sequencing from the dates listed above is stored in bams so the data was converted to fastq outside of the pipeline and stored in the data folder

## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential 1, 2, ... numbering of the gene to a "reference" numbering scheme that represents the standard naming of sites for this gene.
Also assigns each site to a region (domain) of the protein.
So must have columns *sequential_site*, *reference_site*, and *region*.

## Mutation-type classification
[data/mutation_design_classification.csv](data/mutation_design_classification.csv) classifies mutations into the different categories of designed mutations.
Should have columns *sequential_site*, *amino_acid*, and *mutation_type*.

## Neutralization standard barcodes
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) barcodes for the neutralization standards.
Must have columns *barcode* and *name*, giving the barcode and name of this neutralization standard set.

## Barcode runs
[barcode_runs.csv](barcode_runs.csv) must contain the following columns (you can optionally include more):

 - `sample`: sample name, must be unique among barcode runs. Sample name must begin with `<library>-<YYMMDD>` where `<library>` is the library and `<YYMMDD>` is the date. It is recommended (but not enforced) that the full format be `<library>-<YYMMDD>-<description>-<replicate>` where `<description>` is a string description with underscores but no dashes, and `<replicate>` is a number.
 - `library`: name of library, must match a library in the barcode-variant table
 - `date`: date of sequencing, specified in a format parseable to a date by `pandas`.
 - `fastq_R1`: path to one more FASTQ R1 sequencing files, multiple files should be semicolon-delimited

## Configuration for analyzing functional effects of mutations
[func_effects_config.yml](func_effects_config.yml) has the configuration for analyzing functional effects of mutations.
The format is explained within the file.

## Configuration for analyzing antibody escape
[antibody_escape_config.yml](antibody_escape_config.yml) has the configuration for analyzing effects of mutations on escape from antibodies or sera.
The format is explained within the file.