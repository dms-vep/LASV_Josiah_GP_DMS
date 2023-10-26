# Description:
# Python script that adds validation
# sequences and reduced sequence list

# Author:
# Caleb Carr

# Imports
import os
import pandas as pd
from Bio import SeqIO

# Functions
def add_sequences(input_sequences, input_validation_sequence_names, all_sequences):

    print(input_validation_sequence_names)
    
    # Initialize list of curr seqeuences
    list_of_sequences = []
        
    # Get current strains in reduced file
    for curr_fasta in SeqIO.parse(input_sequences, "fasta"):
        curr_name = str(curr_fasta.description)
        list_of_sequences.append(curr_name)

    # Initialize dict for validation sequences
    validation_seqs = dict()

    # Get validation sequences from all sequences
    for curr_fasta in SeqIO.parse(all_sequences, "fasta"):
        curr_name = str(curr_fasta.description)
        if curr_name.split(" ")[0][:-2] in input_validation_sequence_names:
            validation_seqs[curr_name] = str(curr_fasta.seq)

    with open(input_sequences, "a+") as output_sequences:
        for strain, sequence in validation_seqs.items():
            # Write current fasta sequence to output file
            # if strain is not present in file
            if strain not in list_of_sequences:
                output_sequences.write(f">{strain}\n")
                output_sequences.write(f"{sequence}\n")

    # Close files
    output_sequences.close()


def main():
    """
    Main method
    """

    # Input files
    input_sequences = str(snakemake.input.fasta_sequences)
    input_validation_sequence_names = snakemake.params.validation_sequence_names
    all_sequences = str(snakemake.input.all_sequences)
    # Output files
    output_log = str(snakemake.output.output_log)

    log_text = ""
    add_sequences(input_sequences, input_validation_sequence_names, all_sequences)
    log_text += "Validation data added!"

    # Write output log file
    output_log_file = open(output_log, "w")
    output_log_file.write(f"{log_text}\n")

    # Close files
    output_log_file.close()


if __name__ == "__main__":
    main()