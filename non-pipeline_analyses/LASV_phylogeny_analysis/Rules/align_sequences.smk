# Description:
# Aligns fasta sequences

# Author:
# Caleb Carr

rule reduce_similar_protein_sequences:
    """
    This rule reduces the number of sequences 
    by removing similar protein sequences.
    """
    input:
        fasta_sequences = config["Protein_sequences"],
    output:
        fasta_reduced_sequences = config["Protein_reduced_sequences"],
    shell:
        # -i    input input filename in fasta format, required
        # -o    output filename, required
        # -c    sequence identity threshold, default 0.9
        #       this is the default cd-hit's "global sequence identity"
        #       calculated as:
        #        number of identical amino acids in alignment
        #        divided by the full length of the shorter sequence
        # -n    word_length, default 5, see user's guide for choosing it
        #        -n 5 for thresholds 0.7 ~ 1.0
        #        -n 4 for thresholds 0.6 ~ 0.7
        #        -n 3 for thresholds 0.5 ~ 0.6
        #        -n 2 for thresholds 0.4 ~ 0.5
        # -s    length difference cutoff, default 0.0
        #       if set to 0.9, the shorter sequences need to be
        #       at least 90% length of the representative of the cluster
        "cd-hit -i {input.fasta_sequences} -o {output.fasta_reduced_sequences} -c 0.975 -n 5 -s 0.99"

rule add_validation_sequences:
    """
    This rule adds the validation seqeunces 
    for antibody escape to list of sequences.
    """
    input:
        fasta_sequences = config["Protein_reduced_sequences"],
        all_sequences = config["Protein_sequences"],
    params:
        validation_sequence_names = config["Validation_sequence_names"],
    output:
        output_log = config["Validation_addition_log"],
    script:
        "../Scripts/add_validation_sequences.py"


rule align_reduced_protein_sequences:
    """
    This rule aligns reduced number of GPC protein sequences using MAFFT
    """
    input:
        fasta_sequences = config["Protein_reduced_sequences"],
        output_log = config["Validation_addition_log"],
    output:
        alignment = config["GPC_protein_reduced_alignment_temp"],
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        # --retree      Guide tree is built number times in the progressive 
        #               stage. Valid with 6mer distance. Default: 2
        # --maxiterate  Number cycles of iterative refinement are performed. Default: 0
        # --quiet       Do not report progress. Default: off
        "mafft --auto {input.fasta_sequences} > {output.alignment}"


rule align_segment_sequences:
    """
    This rule aligns all S segment sequences using MAFFT
    """
    input:
        fasta_sequences = config["Nucleotide_sequences"],
        output_log = config["Outgroup_addition_log"],
    output:
        alignment = config["S_segment_alignment"],
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        # --retree      Guide tree is built number times in the progressive 
        #               stage. Valid with 6mer distance. Default: 2
        # --maxiterate  Number cycles of iterative refinement are performed. Default: 0
        # --quiet       Do not report progress. Default: off
        "mafft --auto {input.fasta_sequences} > {output.alignment}"


rule align_protein_sequences:
    """
    This rule aligns all GPC protein sequences using MAFFT
    """
    input:
        fasta_sequences = config["Protein_sequences"],
    output:
        alignment = config["GPC_protein_alignment_temp"],
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        # --retree      Guide tree is built number times in the progressive 
        #               stage. Valid with 6mer distance. Default: 2
        # --maxiterate  Number cycles of iterative refinement are performed. Default: 0
        # --quiet       Do not report progress. Default: off
        "mafft --auto {input.fasta_sequences} > {output.alignment}"


rule create_codon_alignment:
    """
    This rule creates a codon aligmnet from the
    codon sequences and protein alignment.
    """
    input: 
        protein_alignment = config["GPC_protein_alignment_temp"],
        codon_sequences = config["Codon_sequences"],
    output:
        config["GPC_codon_alignment"]
    script:
        "../Scripts/create_codon_alignment.py"


rule create_reduced_codon_alignment:
    """
    This rule creates a codon aligmnet from the
    codon sequences and protein alignment.
    """
    input: 
        protein_alignment = config["GPC_protein_reduced_alignment_temp"],
        codon_sequences = config["Codon_sequences"],
    output:
        config["GPC_codon_reduced_alignment"]
    script:
        "../Scripts/create_codon_alignment.py"


rule edit_fasta_headers:
    """
    This rule edits the protein alignment fasta
    headers to remove the appended info from EMBOSS
    getorf.
    """
    input: 
        protein_alignment = config["GPC_protein_alignment_temp"],
    output:
        protein_alignment = config["GPC_protein_alignment"],
    script:
        "../Scripts/edit_fasta_headers.py"


rule edit_reduced_fasta_headers:
    """
    This rule edits the protein alignment fasta
    headers to remove the appended info from EMBOSS
    getorf.
    """
    input: 
        protein_alignment = config["GPC_protein_reduced_alignment_temp"],
    output:
        protein_alignment = config["GPC_protein_reduced_alignment"],
    script:
        "../Scripts/edit_fasta_headers.py"
