# Description:
# Python script that removes the extra
# info appended to fasta headers from EMBOSS
# getorf

# Author:
# Caleb Carr

# Imports
import os
from Bio import SeqIO

# Functions
def edit_fasta_headers(input_protein_alignment, output_protein_alignment):

    # Iterate through headers and remove added info
    with open(output_protein_alignment, "w") as new_protein_alignment:
        for curr_fasta in SeqIO.parse(input_protein_alignment, "fasta"):
            new_name = str(curr_fasta.description).split(" ")[0][:-2]
            curr_fasta.id = new_name
            curr_fasta.description = ""
            SeqIO.write(curr_fasta, new_protein_alignment, "fasta")

    # Close files
    new_protein_alignment.close()

def main():
    """
    Main method
    """

    # Input files
    input_protein_alignment = str(snakemake.input.protein_alignment)

    # Output files
    output_protein_alignment = str(snakemake.output.protein_alignment)

    edit_fasta_headers(input_protein_alignment, output_protein_alignment)


if __name__ == "__main__":
    main()