# Description:
# Python script that removes any sequences 
# that have duplicates in the alignment fasta

# Author:
# Caleb Carr

# Imports
import os
from Bio import SeqIO

# Functions
def remove_duplicate_sequences(input_alignment_file, output_alignment_file, outgroup, log_output_file):

    # Initialize set of seen sequences
    seen = set()
    total = 0
    retained = 0

    # Iterate through headers and remove added info
    with open(output_alignment_file, "w") as new_alignment:
        for curr_fasta in SeqIO.parse(input_alignment_file, "fasta"):
            total += 1
            if curr_fasta.description == outgroup or curr_fasta.seq not in seen:
                seen.add(curr_fasta.seq)
                retained += 1
                SeqIO.write(curr_fasta, new_alignment, "fasta")

    # Write output log file
    output_log_file = open(log_output_file, "a")
    output_log_file.write(
        f"{input_alignment_file} processed a total of "
        f"{total} sequences and retained "
        f"{retained} while removing "
        f"{total-retained} duplicates!\n"
    )

    # Close files
    output_log_file.close()
    new_alignment.close()

def main():
    """
    Main method
    """

    # Input files
    input_segment_alignment = str(snakemake.input.segment_alignment)
    input_codon_alignment = str(snakemake.input.codon_alignment)
    input_protein_alignment = str(snakemake.input.protein_alignment)
    # Params
    outgroup = str(snakemake.params.outgroup)

    # Output files
    log_output_file = str(snakemake.output.log)
    output_segment_alignment = str(snakemake.output.segment_alignment)
    output_codon_alignment = str(snakemake.output.codon_alignment)
    output_protein_alignment = str(snakemake.output.protein_alignment)

    remove_duplicate_sequences(input_segment_alignment, output_segment_alignment, outgroup, log_output_file)
    remove_duplicate_sequences(input_codon_alignment, output_codon_alignment, outgroup, log_output_file)
    remove_duplicate_sequences(input_protein_alignment, output_protein_alignment, outgroup, log_output_file)



if __name__ == "__main__":
    main()