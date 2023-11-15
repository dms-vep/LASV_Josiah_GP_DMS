# Description:
# Aligns fasta sequences

# Author:
# Caleb Carr


rule align_reduced_protein_sequences:
    """
    This rule aligns reduced number of protein sequences using MAFFT
    """
    input:
        fasta_sequences = config["Protein_important_sequences"],
    output:
        alignment = config["Reduced_protein_alignment"],
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
    This rule aligns all protein sequences using MAFFT
    """
    input:
        fasta_sequences = config["Protein_sequences"],
    output:
        alignment = config["Protein_alignment"],
    shell:
        # --auto        Automatically selects an appropriate strategy from L-INS-i, 
        #               FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        # --retree      Guide tree is built number times in the progressive 
        #               stage. Valid with 6mer distance. Default: 2
        # --maxiterate  Number cycles of iterative refinement are performed. Default: 0
        # --quiet       Do not report progress. Default: off
        "mafft --auto {input.fasta_sequences} > {output.alignment}"
