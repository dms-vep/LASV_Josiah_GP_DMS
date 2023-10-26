# Description:
# Python script that creates a codon alignment from 
# the protein alignment and codon fasta sequences

# Author:
# Caleb Carr

# Imports
from Bio import SeqIO

# Functions
def create_codon_alignment(protein_alignment_file, codon_fasta_file, codon_alignemnt_file):

    # Initialize dicts for sequences
    protein_alignment_dict = {}
    codon_alignment_dict = {}

    # Get protein alignment sequences
    for curr_seq in SeqIO.parse(protein_alignment_file, "fasta"):
        name = curr_seq.description
        protein_alignment_dict[name] = str(curr_seq.seq)

    # Get codon fasta sequences
    for curr_seq in SeqIO.parse(codon_fasta_file, "fasta"):
        name = curr_seq.description
        codon_alignment_dict[name] = str(curr_seq.seq)

    # Iterate through protein alignments and creating codon alignments
    with open(codon_alignemnt_file, "w") as codon_alignment:
        for header, sequence in protein_alignment_dict.items():
            curr_codon_seq = codon_alignment_dict[header]
            codon_alignment_seq = ""
            codon_index = 0
            for amino_acid in sequence:
                if amino_acid == "-":
                    codon_alignment_seq += "---"
                else:
                    codon_alignment_seq += curr_codon_seq[codon_index:codon_index + 3]
                    codon_index += 3

            codon_alignment.write(f">{header.split(' ')[0][:-2]}\n")
            codon_alignment.write(f"{codon_alignment_seq}\n")

    # Close files
    codon_alignment.close()

def main():
    """
    Main method
    """

    # Input files
    protein_alignment_file = str(snakemake.input.protein_alignment)
    codon_fasta_file = str(snakemake.input.codon_sequences)

    # Output files
    codon_alignemnt_file = str(snakemake.output)

    create_codon_alignment(protein_alignment_file, codon_fasta_file, codon_alignemnt_file)


if __name__ == "__main__":
    main()