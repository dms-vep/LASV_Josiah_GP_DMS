# Non-pipeline analyses
The code and data used for library construction, phylogeny analysis, virus titers, virus neutralization assays, etc. are included in the following subdirectories.

## Repo organization
 
* `DAG1_phylogeny_analysis` directory contains a fully automated snakemake pipeline that analyzes DAG1 vertebrate sequences. 
* `LASV_NGS_analysis` directory contains a fully automated snakemake pipeline that analyzes sequencing data for the LM395 strain. 
* `LASV_phylogeny_analysis` directory contains a fully automated snakemake pipeline that analyzes all available Lassa virus strain sequences. 
* `library_construction` directory contains directories for library design and initial PCR mutagenesis QC (using Sanger sequencing).
* `pdb_antibody_contacts` directory contains a script to calculate distances between GPC sites and antibody contact sites. 
* `plasmid_maps` directory contains all plasmid maps for plasmids used in this study.
* `primers` directory contains all primer sequences for primers used in this study.  
* `titers_and_neuts` directory contains scripts and data for testing virus titers and monoclonal antibody neutralization. The data includes a combination of flow-based and luciferase-based results.
