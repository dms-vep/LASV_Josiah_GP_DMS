# Description:
# Downloads accessions and process/extract sequences and metadata

# Author:
# Caleb Carr


rule process_accession_list:
    """
    This rule extracts lists of matching mRNA and protein accessions
    """
    input:
        config["Ortholog_accession_csv"],
    output:
        mRNA_accession_list = config["mRNA_accession_list"],
        protein_accession_list = config["Protein_accession_list"],
    script:
        "../Scripts/process_accession_lists.py"

rule download_and_process_accessions:
    """
    This rule runs the download_NCBI_sequences.py script which
    downloads all genbank files from a list of accessions and 
    extracts metadata. 
    """
    input:
        mRNA_accession_list = config["mRNA_accession_list"],
        protein_accession_list = config["Protein_accession_list"],
    params:
        genome_size_threshold_lower = config["Genome_size_threshold_lower"],
        genome_size_threshold_upper = config["Genome_size_threshold_upper"],
    output:
        CDS_fasta_sequences = config["CDS_sequences"],
        protein_fasta_sequences = config["Protein_sequences"],
        metadata = config["Metadata"],
    script:
        "../Scripts/download_NCBI_sequences.py"