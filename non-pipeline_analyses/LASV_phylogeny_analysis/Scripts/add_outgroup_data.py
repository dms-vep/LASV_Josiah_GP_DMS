# Description:
# Python script that adds outgroup
# sequences and metadata to existing
# files

# Author:
# Caleb Carr

# Imports
import os
import pandas as pd
from Bio import SeqIO

# Functions
def add_metadata(input_metadata, input_outgroup_metadata, output_metadata):

    # Load metadata as dataframe
    metadata = pd.read_csv(input_metadata, sep="\t")
    outgroup_metadata = pd.read_csv(input_outgroup_metadata, sep="\t")

    # Append outgroup data
    metadata = (
        pd.concat([
            metadata, 
            outgroup_metadata,
            ], ignore_index = True)
    )

    # Write updated metadata to file
    metadata.to_csv(output_metadata, sep="\t", index=False)

def add_sequences(input_sequences, input_outgroup_sequences, output_sequences):

    with open(output_sequences, "w") as output_sequences:
        # Write current sequences to new file
        for curr_fasta in SeqIO.parse(input_sequences, "fasta"):
            strain = str(curr_fasta.id)
            sequence = str(curr_fasta.seq)
            # Write current fasta sequence to output file
            output_sequences.write(f">{strain}\n")
            output_sequences.write(f"{sequence}\n")
        # Write outgroup sequences to new file
        for curr_fasta in SeqIO.parse(input_outgroup_sequences, "fasta"):
            strain = str(curr_fasta.id)
            sequence = str(curr_fasta.seq)
            # Write current fasta sequence to output file
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
    input_metadata = str(snakemake.input.metadata)
    input_outgroup_sequences = str(snakemake.input.outgroup_fasta_sequences)
    input_outgroup_metadata = str(snakemake.input.outgroup_metadata)
    # Output files
    output_sequences = str(snakemake.output.fasta_sequences)
    output_metadata = str(snakemake.output.metadata)

    add_metadata(input_metadata, input_outgroup_metadata, output_metadata)
    add_sequences(input_sequences, input_outgroup_sequences, output_sequences)


if __name__ == "__main__":
    main()