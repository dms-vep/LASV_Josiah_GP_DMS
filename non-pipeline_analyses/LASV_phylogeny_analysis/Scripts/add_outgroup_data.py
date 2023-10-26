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
def add_metadata(input_metadata, input_outgroup_metadata):

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
    metadata.to_csv(input_metadata, sep="\t", index=False)

def add_sequences(input_sequences, input_outgroup_sequences):

    with open(input_sequences, "a+") as output_sequences:
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
    # Params
    add_outgroup = snakemake.params.add_outgroups
    # Output files
    output_log = str(snakemake.output.output_log)

    log_text = ""

    if add_outgroup == True:
        add_metadata(input_metadata, input_outgroup_metadata)
        add_sequences(input_sequences, input_outgroup_sequences)
        log_text += "Outgroup data added!"
    else:
        log_text += "Outgroup data not added!"

    # Write output log file
    output_log_file = open(output_log, "w")
    output_log_file.write(f"{log_text}\n")

    # Close files
    output_log_file.close()


if __name__ == "__main__":
    main()