# Description:
# Processes nucleotide fastas to find largest ORF in 
# forward direction and translates to protein sequences

# Author:
# Caleb Carr


rule get_protein_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    amino acid sequences for GPC segment. 
    """
    input:
        sequences = config["Nucleotide_sequences"],
        outgroup_log = config["Outgroup_addition_log"],
    params:
        config["Min_ORF_threshold"],
    output:
        config["Protein_sequences"],
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 1' flag means the amino acid 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 1 -minsize {params} -reverse No"


rule get_codon_sequences:
    """
    This rule runs EMBOSS getorf to extract 
    codon sequences for GPC segment. 
    """
    input:
        sequences = config["Nucleotide_sequences"],
        outgroup_log = config["Outgroup_addition_log"],
    params:
        config["Min_ORF_threshold"],
    output:
        config["Codon_sequences"],
    shell:
        # The '-sequence' flag signals for the input file
        # while the '-outseq' file signals for the output 
        # file. The '-find 3' flag means the nucleotide 
        # sequences between START and STOP codons are returned.
        # The '-minsize' flag means minimum nucleotide size of 
        # ORF to report. The '-reverse' flag signals if ORFs on 
        # the reverse strand should be found as well. 
        "getorf -sequence {input.sequences} -outseq {output} -find 3 -minsize {params} -reverse No"


rule check_number_ORFs_found:
    """
    This rule checks both the amino acid
    and codon fastas to verify there is a 
    one-to-one mapping from S segment nucleotide
    seqeunce to protein and codon sequences. 
    """
    input:
        nucleotide = config["Nucleotide_sequences"],
        protein = config["Protein_sequences"],
        codon = config["Codon_sequences"],
    output:
        config["Sequence_verification_log"]
    script:
        "../Scripts/extracted_sequence_verification.py"