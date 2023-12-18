# Description:
# Downloads accessions and process/extract sequences and metadata

# Author:
# Caleb Carr


rule download_and_process_accessions:
    """
    This rule runs the download_NCBI_sequences.py script which
    downloads all genbank files from a list of accessions and 
    extracts metadata. 
    """
    input:
        accession_list = config["Accession_list"],
        reference_genome = config["Reference_genome"],
    params:
        genome_size_threshold_lower = config["Genome_size_threshold_lower"],
        genome_size_threshold_upper = config["Genome_size_threshold_upper"],
        desired_segment = config["Desired_segment"],
        temp_fasta_file = config["Temp_fasta_file"],
        temp_alignment_file = config["Temp_alignment_file"],
        max_frac_N = config["max_frac_N"],
        accesstions_to_exclude = config["Accessions_to_exclude"],
    output:
        fasta_sequences = config["Nucleotide_sequences"],
        metadata = config["Metadata"],
    script:
        "../Scripts/download_NCBI_sequences.py"


rule add_outgroup_references:
    """
    This rule adds the outgroup strains to the
    fasta files and metadata files.
    """
    input:
        fasta_sequences = config["Nucleotide_sequences"],
        metadata = config["Metadata"],
        outgroup_fasta_sequences = config["Outgroup_reference_genomes"],
        outgroup_metadata = config["Outgroup_reference_metadata"],
    output:
        fasta_sequences = config["Nucleotide_sequences_with_outgroup"],
        metadata = config["Metadata_with_outgroup"],
    script:
        "../Scripts/add_outgroup_data.py"
