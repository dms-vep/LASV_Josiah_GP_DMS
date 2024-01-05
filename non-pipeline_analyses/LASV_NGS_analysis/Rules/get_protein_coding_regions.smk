
# Description:
# Extract largest (assumed to be the main formed protein) 
# protein coding region from each gene.

# Author:
# Caleb Carr


rule python_create_consensus_fasta_sequence:
    """ This rule creates a consensus sequence for each sample 
        in fasta format based on the mpileup summary. """
    input:
        "Results/Sample_Data/{accession}/{accession}_mpileup_summary.csv",
    output:
        "Results/Sample_Data/{accession}/{accession}.fasta",
    params:
        gene_list = config['Gene_list'],
    conda:
        "../Conda_Envs/consensus_protein_env.yml"
    script:
        # Python script to get consensus sequence from
        # summary file
        "../Scripts/get_consensus_sequence_fasta.py"


rule emboss_get_all_possible_orfs:
    """ This rule gets all possible orfs for each gene segement 
        based on the newly created fasta sequence using EMBOSS getorf. """
    input:
        "Results/Sample_Data/{accession}/{accession}.fasta",
    output:
        temp("Results/Sample_Data/{accession}/{accession}_all_orfs.fasta"),
    conda:
        "../Conda_Envs/consensus_protein_env.yml"
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 3' flag means the nucleotide 
        # sequences between START and STOP codons are returned.
        "getorf -sequence {input} -outseq {output} -find 3"


rule python_get_most_likely_protein_coding_sequence:
    """ This rule finds the largest segments which is assumed to 
        correspond to the major protein products and output the 
        result as a fasta file with just the protein coding segments 
        by parsing the file with all orfs. """
    input:
        "Results/Sample_Data/{accession}/{accession}_all_orfs.fasta",
    output:
        sequences = "Results/Sample_Data/{accession}/{accession}_protein.fasta",
        positions = "Results/Sample_Data/{accession}/{accession}_protein_indices.csv",
    conda:
        "../Conda_Envs/consensus_protein_env.yml"
    script:
        # Python script to get largest orfs for each gene
        "../Scripts/get_largest_orfs.py"